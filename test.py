import os
import sys
import subprocess
import sqlite3
import pickle
import copy

import matplotlib.colors as mcolors
from HiPRGen.network_loader import NetworkLoader
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge
from monty.serialization import loadfn, dumpfn
from HiPRGen.species_filter import species_filter, add_electron_species
from HiPRGen.bucketing import bucket
from HiPRGen.report_generator import ReportGenerator
from HiPRGen.initial_state import insert_initial_state
from HiPRGen.constants import ROOM_TEMP, Terminal
from HiPRGen.reaction_filter_payloads import DispatcherPayload, WorkerPayload

from HiPRGen.species_questions import (
    mg_species_decision_tree,
    li_species_decision_tree,
    nonmetal_species_decision_tree,
    euvl_species_decision_tree,
    positive_penalty,
    species_default_true,
)

from HiPRGen.reaction_questions import (
    default_reaction_decision_tree,
    default_logging_decision_tree,
    co2_reaction_decision_tree,
    euvl_phase1_reaction_decision_tree,
    euvl_phase1_logging_decision_tree,
    euvl_phase2_reaction_decision_tree,
    euvl_phase2_logging_decision_tree
)

from HiPRGen.mc_analysis import (
    reaction_tally_report,
    species_report,
    Pathfinding,
    SimulationReplayer,
    pathway_report,
    sink_report,
    consumption_report,
    redox_report,
    final_state_report,
)

# Since HiPRGen uses an end-to-end testing approach rather than testing
# each individual function, we have decided to use the tests as
# documentation, by explaining every single line through the first test.


# The first thing you need to consider when using HiPRGen is how many
# worker threads do you want to run. HiPRGen can be run with a single
# thread or thousands distrubuted across several nodes. For reaction
# networks with between ~5000 and ~10000 species, we have found that the
# optimal number of worker threads is between 1000 and 2000. If you try
# and use more than that, the worker threads are going to spend lots of
# time waiting for the dispatcher to get through all of the reactions it
# is being sent, which slows everything down. Fixing this would require
# a more complex distrubuted system, but it hasn't been an issue even
# for very large reaction networks.
if len(sys.argv) != 2:
    print("usage: python test.py number_of_threads")
    quit()


number_of_threads = sys.argv[1]


class bcolors:
    PASS = "\u001b[32;1m"
    FAIL = "\u001b[31;1m"
    ENDC = "\u001b[0m"


# HiPRGen is organized as a pipeline, where all the relevent data is
# stored in a sqlite database between phases. For this reason, during
# a single run of the full pipeline, it makes sense to store all the
# relevent files in a single directory. We have two test sets, a lithium
# set and a magnesium set. Since the lithium test set is older, we shall
# document that instead of the mg test set.

if os.path.isdir("./scratch"):
    subprocess.run(["rm", "-r", "./scratch"])

subprocess.run(["mkdir", "./scratch"])


def li_test():

    # folder is the where we store all our intermediate databases
    folder = "./scratch/li_test"
    subprocess.run(["mkdir", folder])

    # The initial input to the pipeline is a list of LIBE or MADEIRA
    # dataset entries. We provide two examples in the data foloder.
    mol_json = "./data/ronald_LIBE.json"
    database_entries = loadfn(mol_json)
    # The first step of the HiPRGen pipeline is passing the input molecules
    # through the species decision tree to discard molecules.
    species_decision_tree = li_species_decision_tree

    params = {"temperature": ROOM_TEMP, "electron_free_energy": -1.4}

    # There is one non-local part of species filtering: we consider two
    # molecules to be equivalent if they have the same total charge,
    # composition, and covalent bonds, even if they have different metal
    # coordination, and we choose one such molecule in each "coordimer"
    # class using the coodimer weight function.  Since most of our logging
    # later on is defined in terms of a fixed molecule set, logging for
    # the species filtering phase is messy, so ignore the species_report
    # argument for now. The second argument is where we store a pickle of
    # the filtered molecule entries for use in later phases.

    mol_entries = species_filter(
        database_entries,
        mol_entries_pickle_location=folder + "/mol_entries.pickle",
        species_report=folder + "/unfiltered_species_report.tex",
        species_decision_tree=species_decision_tree,
        coordimer_weight=lambda mol: (mol.penalty, mol.solvation_correction + mol.get_free_energy(params["temperature"])),
    )

    # Once we have generated our molecule list, we generate the bucket database
    # which is how we break up the reaction filtering amongst all avaliable workers.
    # It gets stored in the buckets.sqlite database.
    bucket(mol_entries, folder + "/buckets.sqlite")

    # Reaction filtering is paralellized using MPI, so we need to spawn
    # an MPI instance to run it. This is why we can't just start
    # reaction filtering by calling a python function. We pass the
    # reaction decision tree, the logging decision tree, and the electron
    # free energy as strings across this barrier. Every possible
    # reaction gets passed through both the reaction decision tree and
    # the logging decision tree. If a reaction passes the reaction
    # decision tree, it gets written to the network. If a reaction
    # passes the logging decision tree, it gets logged to the reaction
    # report along with what happened to it in reaction_decision_tree.

    # The reaction decision trees are constructed in
    # HiPRGen.reaction_questions

    dispatcher_payload = DispatcherPayload(
        folder + "/buckets.sqlite",
        folder + "/rn.sqlite",
        folder + "/reaction_report.tex",
    )

    worker_payload = WorkerPayload(
        folder + "/buckets.sqlite",
        default_reaction_decision_tree,
        params,
        default_logging_decision_tree,
    )

    # The dispatcher and worker payloads are passed through the MPI barrier
    # as JSON blobs dispatcher_payload and worker_payload
    dumpfn(dispatcher_payload, folder + "/dispatcher_payload.json")
    dumpfn(worker_payload, folder + "/worker_payload.json")

    subprocess.run(
        [
            "mpirun",
            "--use-hwthread-cpus",
            "-n",
            number_of_threads,
            "python",
            "run_network_generation.py",
            folder + "/mol_entries.pickle",
            folder + "/dispatcher_payload.json",
            folder + "/worker_payload.json",
        ]
    )

    # After we have generated the mol_entries, we refer to molecules by
    # their index. The function find_mol_entry_from_xyz_and_charge can
    # help find the indices of specific species to be used in the initial
    # condition for propagating trajectories and/or trajectory analysis.
    Li_plus_id = find_mol_entry_from_xyz_and_charge(
        mol_entries, "./xyz_files/Li.xyz", 1
    )

    EC_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/EC.xyz", 0)

    LEDC_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/LEDC.xyz", 0)

    # After generating a reaction network, it is stored in rn.sqlite. We
    # use Monte Carlo simulation to interrogate the network, and for that
    # we need to define an initial condition.
    initial_state = {Li_plus_id: 30, EC_id: 30}

    # The initial state and the trajectories (after simulation) are stored in
    # a seperate database from the network, here called initial_state.sqlite.
    # This facilitates running multiple independent simulations of the same
    # network with different initial conditions at the same time, if desired.
    insert_initial_state(initial_state, mol_entries, folder + "/initial_state.sqlite")

    # GMC is a high performance reaction network Monte Carlo simulator using the
    # Gillespie algorithm: https://github.com/BlauGroup/RNMC. Here we run 1000
    # trajectories each of 200 steps.
    subprocess.run(
        [
            "GMC",
            "--reaction_database=" + folder + "/rn.sqlite",
            "--initial_state_database=" + folder + "/initial_state.sqlite",
            "--number_of_simulations=1000",
            "--base_seed=1000",
            "--thread_count=" + number_of_threads,
            "--step_cutoff=200",
        ]
    )

    # The network loader builds a python object around a reaction network
    # and the molecules to make it easier to use them.
    network_loader = NetworkLoader(
        folder + "/rn.sqlite",
        folder + "/mol_entries.pickle",
        folder + "/initial_state.sqlite",
    )

    network_loader.load_initial_state_and_trajectories()

    # HiPRGen has analysis tools to understand what happened in our simulation.
    # The output files are written into the same folder in which the reaction
    # network is stored.

    # This report is empty, but we use it to generate the molecule pictures.
    # This is an expensive operation, so we only want do do it once.
    report_generator = ReportGenerator(
        network_loader.mol_entries, folder + "/dummy.tex", rebuild_mol_pictures=True
    )

    # The tally report shows reactions sorted by the number of times fired.
    reaction_tally_report(network_loader, folder + "/reaction_tally.tex")
    # Run `pdflatex reaction_tally.tex` in `scratch/li_test` to generate
    # the tally report PDF.

    # The species report shows every specie in the network and their IDs.
    species_report(network_loader, folder + "/species_report.tex")
    # Run `pdflatex species_report.tex` in `scratch/li_test` to generate
    # the species report PDF.

    # Pathfinding is a central goal of HiPRGen / GMC. See mc_analysis.py for
    # further documentation of the Pathfinding class.
    pathfinding = Pathfinding(network_loader)

    # The pathway report shows all the ways that a target species was
    # produced in the simulation trajectories, where each simulation only
    # contributes the shortest path responsible for the first formation
    # of the target species to the report. The report can be sorted by
    # pathway frequency, but instead here we sort by pathway cost. Note
    # that the test network has ~5000 reactions while production networks
    # have between 50-100 million reactions.
    pathway_report(
        pathfinding, LEDC_id, folder + "/LEDC_pathways.tex", sort_by_frequency=False
    )
    # Run `pdflatex LEDC_pathways.tex` in `scratch/li_test` to generate
    # the LEDC pathway report PDF.

    # The simulation replayer sweeps through all trajectories in order
    # to extract additional information that is used for consumption
    # reports and sink reports.
    simulation_replayer = SimulationReplayer(network_loader)

    # The consumption report shows reactions which consumed a target
    # species, sorted by the number of times the reaction fired.
    consumption_report(
        simulation_replayer, LEDC_id, folder + "/LEDC_consumption_report.tex"
    )
    # Run `pdflatex LEDC_consumption_report.tex` in `scratch/li_test`
    # to generate the LEDC consumption report PDF.

    # The sink report shows species which have a production to
    # consumption ratio of greater than 3/2 and which have an expected
    # value above 0.1. These are two of the three heuristic criteria
    # that we use to identify network products. The third criteria is
    # that each network product must have a shortest path with cost
    # less than 10. This can be checked by generating pathway reports
    # to each species shown in the sink report. For the curious reader,
    # we note that generating pathway reports to the six species in the
    # sink report will show that only Li2CO3, C2H4, LiEDC-, and DLEMC
    # have sufficiently low-cost paths to pass the third criteria and
    # thus to be considered products of the test network used here.
    sink_report(simulation_replayer, folder + "/sink_report.tex")
    # Run `pdflatex sink_report.tex` in `scratch/li_test` to generate
    # the sink report PDF.

    tests_passed = True
    print("Number of species:", network_loader.number_of_species)
    if network_loader.number_of_species == 190:
        print(bcolors.PASS + "li_test: correct number of species" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "li_test: correct number of species" + bcolors.ENDC)
        tests_passed = False

    print("Number of reactions:", network_loader.number_of_reactions)
    if network_loader.number_of_reactions == 4649:
        print(bcolors.PASS + "li_test: correct number of reactions" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "li_test: correct number of reactions" + bcolors.ENDC)
        tests_passed = False

    return tests_passed


def mg_test():

    folder = "./scratch/mg_test"
    subprocess.run(["mkdir", folder])

    mol_json = "./data/sam_G2.json"
    species_decision_tree = mg_species_decision_tree

    database_entries = loadfn(mol_json)
    params = {"temperature": ROOM_TEMP, "electron_free_energy": -2.06}

    mol_entries = species_filter(
        database_entries,
        folder + "/mol_entries.pickle",
        folder + "/unfiltered_species_report.tex",
        species_decision_tree,
        coordimer_weight=lambda mol: (mol.penalty, mol.solvation_correction + mol.get_free_energy(params["temperature"])),
    )

    bucket(mol_entries, folder + "/buckets.sqlite")

    dispatcher_payload = DispatcherPayload(
        folder + "/buckets.sqlite",
        folder + "/rn.sqlite",
        folder + "/reaction_report.tex",
    )

    worker_payload = WorkerPayload(
        folder + "/buckets.sqlite",
        default_reaction_decision_tree,
        params,
        default_logging_decision_tree,
    )

    dumpfn(dispatcher_payload, folder + "/dispatcher_payload.json")
    dumpfn(worker_payload, folder + "/worker_payload.json")

    subprocess.run(
        [
            "mpiexec",
            "--use-hwthread-cpus",
            "-n",
            number_of_threads,
            "python",
            "run_network_generation.py",
            folder + "/mol_entries.pickle",
            folder + "/dispatcher_payload.json",
            folder + "/worker_payload.json",
        ]
    )

    mg_g2_plus_plus_id = find_mol_entry_from_xyz_and_charge(
        mol_entries, "./xyz_files/mgg2.xyz", 2
    )

    c2h4_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/c2h4.xyz", 0)

    c2h6_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/c2h6.xyz", 0)

    initial_state = {33: 30, 81: 30}

    insert_initial_state(initial_state, mol_entries, folder + "/initial_state.sqlite")

    subprocess.run(
        [
            "GMC",
            "--reaction_database=" + folder + "/rn.sqlite",
            "--initial_state_database=" + folder + "/initial_state.sqlite",
            "--number_of_simulations=1000",
            "--base_seed=1000",
            "--thread_count=" + number_of_threads,
            "--step_cutoff=200",
        ]
    )

    network_loader = NetworkLoader(
        folder + "/rn.sqlite",
        folder + "/mol_entries.pickle",
        folder + "/initial_state.sqlite",
    )

    network_loader.load_initial_state_and_trajectories()

    report_generator = ReportGenerator(
        network_loader.mol_entries, folder + "/dummy.tex", rebuild_mol_pictures=True
    )

    reaction_tally_report(network_loader, folder + "/reaction_tally.tex")

    pathfinding = Pathfinding(network_loader)

    pathway_report(
        pathfinding, c2h6_id, folder + "/C2H6_pathways.tex", sort_by_frequency=False
    )

    pathway_report(
        pathfinding, c2h4_id, folder + "/C2H4_pathways.tex", sort_by_frequency=False
    )

    species_report(network_loader, folder + "/species_report.tex")

    tests_passed = True
    print("Number of species:", network_loader.number_of_species)
    if network_loader.number_of_species == 83:
        print(bcolors.PASS + "mg_test: correct number of species" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "mg_test: correct number of species" + bcolors.ENDC)
        tests_passed = False

    print("Number of reactions:", network_loader.number_of_reactions)
    if network_loader.number_of_reactions == 766:
        print(bcolors.PASS + "mg_test: correct number of reactions" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "mg_test: correct number of reactions" + bcolors.ENDC)
        tests_passed = False

    return tests_passed


# def flicho_test():

#     folder = "./scratch/flicho_test"
#     subprocess.run(["mkdir", folder])

#     mol_json = "./data/flicho_test.json"
#     database_entries = loadfn(mol_json)
#     species_decision_tree = li_species_decision_tree

#     mol_entries = species_filter(
#         database_entries,
#         mol_entries_pickle_location=folder + "/mol_entries.pickle",
#         species_report=folder + "/unfiltered_species_report.tex",
#         species_decision_tree=species_decision_tree,
#         coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy),
#     )

#     bucket(mol_entries, folder + "/buckets.sqlite")

#     params = {"temperature": ROOM_TEMP, "electron_free_energy": -1.4}

#     dispatcher_payload = DispatcherPayload(
#         folder + "/buckets.sqlite",
#         folder + "/rn.sqlite",
#         folder + "/reaction_report.tex",
#     )

#     worker_payload = WorkerPayload(
#         folder + "/buckets.sqlite",
#         default_reaction_decision_tree,
#         params,
#         Terminal.DISCARD,
#     )

#     dumpfn(dispatcher_payload, folder + "/dispatcher_payload.json")
#     dumpfn(worker_payload, folder + "/worker_payload.json")

#     subprocess.run(
#         [
#             "mpirun",
#             "--use-hwthread-cpus",
#             "-n",
#             number_of_threads,
#             "python",
#             "run_network_generation.py",
#             folder + "/mol_entries.pickle",
#             folder + "/dispatcher_payload.json",
#             folder + "/worker_payload.json",
#         ]
#     )

#     Li_plus_id = find_mol_entry_from_xyz_and_charge(
#         mol_entries, "./xyz_files/Li.xyz", 1
#     )

#     EC_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/EC.xyz", 0)

#     initial_state = {Li_plus_id: 30, EC_id: 30}

#     insert_initial_state(initial_state, mol_entries, folder + "/initial_state.sqlite")

#     subprocess.run(
#         [
#             "GMC",
#             "--reaction_database=" + folder + "/rn.sqlite",
#             "--initial_state_database=" + folder + "/initial_state.sqlite",
#             "--number_of_simulations=1000",
#             "--base_seed=1000",
#             "--thread_count=" + number_of_threads,
#             "--step_cutoff=200",
#         ]
#     )

#     network_loader = NetworkLoader(
#         folder + "/rn.sqlite",
#         folder + "/mol_entries.pickle",
#         folder + "/initial_state.sqlite",
#     )

#     network_loader.load_initial_state_and_trajectories()

#     report_generator = ReportGenerator(
#         network_loader.mol_entries, folder + "/dummy.tex", rebuild_mol_pictures=True
#     )

#     coordination_report(network_loader, folder + "/coodination_report.tex", "Li1", 1)

#     decoordination_report(
#         network_loader, folder + "/decoodination_report.tex", "Li1", 1
#     )

#     return tests_passed


# def co2_test():

#     folder = "./scratch/co2_test"
#     subprocess.run(["mkdir", folder])

#     mol_json = "./data/co2_summary_docs.json"
#     database_entries = loadfn(mol_json)

#     species_decision_tree = nonmetal_species_decision_tree

#     mol_entries = species_filter(
#         database_entries,
#         mol_entries_pickle_location=folder + "/mol_entries.pickle",
#         species_report=folder + "/unfiltered_species_report.tex",
#         species_decision_tree=species_decision_tree,
#         coordimer_weight=lambda mol: (mol.free_energy),
#         species_logging_decision_tree=species_decision_tree,
#         generate_unfiltered_mol_pictures=True,
#     )

#     bucket(mol_entries, folder + "/buckets.sqlite")

#     params = {"temperature": ROOM_TEMP, "electron_free_energy": -4.04}

#     dispatcher_payload = DispatcherPayload(
#         folder + "/buckets.sqlite",
#         folder + "/rn.sqlite",
#         folder + "/reaction_report.tex",
#     )

#     worker_payload = WorkerPayload(
#         folder + "/buckets.sqlite",
#         co2_reaction_decision_tree,
#         params,
#         co2_reaction_decision_tree,
#     )

#     dumpfn(dispatcher_payload, folder + "/dispatcher_payload.json")
#     dumpfn(worker_payload, folder + "/worker_payload.json")

#     subprocess.run(
#         [
#             "mpirun",
#             "--use-hwthread-cpus",
#             "-n",
#             number_of_threads,
#             "python",
#             "run_network_generation.py",
#             folder + "/mol_entries.pickle",
#             folder + "/dispatcher_payload.json",
#             folder + "/worker_payload.json",
#         ]
#     )

    return


def euvl_phase1_test():

    folder = "./scratch/euvl_phase1_test"
    subprocess.run(["mkdir", folder])

    # mol_json = "./data/mini_euvl_test_set.json"
    mol_json = "./data/euvl_test_set.json"
    database_entries = loadfn(mol_json)

    species_decision_tree = euvl_species_decision_tree

    params = {
        "temperature": ROOM_TEMP,
        "electron_free_energy": 0.0,
    }

    mol_entries = species_filter(
        database_entries,
        mol_entries_pickle_location=folder + "/mol_entries.pickle",
        species_report=folder + "/unfiltered_species_report.tex",
        species_decision_tree=species_decision_tree,
        coordimer_weight=lambda mol: (mol.get_free_energy(params["temperature"])),
        species_logging_decision_tree=species_decision_tree,
        generate_unfiltered_mol_pictures=False,
    )

    print(len(mol_entries), "initial mol entries")

    bucket(mol_entries, folder + "/buckets.sqlite")

    params["electron_species"] = len(mol_entries)

    mol_entries = add_electron_species(
        mol_entries,
        mol_entries_pickle_location=folder + "/mol_entries.pickle",
        electron_free_energy=params["electron_free_energy"],
    )

    print(len(mol_entries), "final mol entries")

    dispatcher_payload = DispatcherPayload(
        folder + "/buckets.sqlite",
        folder + "/rn.sqlite",
        folder + "/reaction_report.tex",
    )

    worker_payload = WorkerPayload(
        folder + "/buckets.sqlite",
        euvl_phase1_reaction_decision_tree,
        params,
        euvl_phase1_logging_decision_tree
    )

    dumpfn(dispatcher_payload, folder + "/dispatcher_payload.json")
    dumpfn(worker_payload, folder + "/worker_payload.json")

    subprocess.run(
        [
            "mpirun",
            "--use-hwthread-cpus",
            "-n",
            number_of_threads,
            "python",
            "run_network_generation.py",
            folder + "/mol_entries.pickle",
            folder + "/dispatcher_payload.json",
            folder + "/worker_payload.json",
        ]
    )

    tps_plus1_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/tps.xyz", 1)
    phs_0_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/phs.xyz", 0)
    tba_0_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/tba.xyz", 0)
    nf_minus1_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/nf.xyz", -1)

    initial_state = {tps_plus1_id: 20, phs_0_id: 24, tba_0_id: 36, nf_minus1_id: 20}

    insert_initial_state(initial_state, mol_entries, folder + "/initial_state.sqlite")

    # subprocess.run(
    #     [
    #         "GMC",
    #         "--reaction_database=" + folder + "/rn.sqlite",
    #         "--initial_state_database=" + folder + "/initial_state.sqlite",
    #         "--number_of_simulations=1000",
    #         "--base_seed=1000",
    #         "--thread_count=" + number_of_threads,
    #         "--step_cutoff=200",
    #         "--energy_budget=92",
    #     ]
    # )

    network_loader = NetworkLoader(
        folder + "/rn.sqlite",
        folder + "/mol_entries.pickle",
        folder + "/initial_state.sqlite",
    )

    # network_loader.load_initial_state_and_trajectories()

    report_generator = ReportGenerator(
        network_loader.mol_entries, folder + "/dummy.tex", rebuild_mol_pictures=True
    )

    # reaction_tally_report(network_loader, folder + "/reaction_tally.tex", cutoff=10)
    # species_report(network_loader, folder + "/species_report.tex")
    # simulation_replayer = SimulationReplayer(network_loader)
    # final_state_report(simulation_replayer, folder + "/final_state_report.tex")

    tests_passed = True
    print("Number of species:", network_loader.number_of_species)
    if network_loader.number_of_species == 104:
        print(bcolors.PASS + "euvl_phase_1_test: correct number of species" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "euvl_phase_1_test: correct number of species" + bcolors.ENDC)
        tests_passed = False

    print("Number of reactions:", network_loader.number_of_reactions)
    if network_loader.number_of_reactions == 560:
        print(bcolors.PASS + "euvl_phase_1_test: correct number of reactions" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "euvl_phase_1_test: correct number of reactions" + bcolors.ENDC)
        tests_passed = False

    return tests_passed


def euvl_phase2_test():

    phase1_folder = "./euvl_phase1_test"
    folder = "./scratch/euvl_phase2_test"
    subprocess.run(["mkdir", folder])

    mol_json = "./data/euvl_test_set.json"
    database_entries = loadfn(mol_json)

    species_decision_tree = euvl_species_decision_tree

    params = {
        "temperature": ROOM_TEMP+200.0,
        "electron_free_energy": 0.0,
    }

    mol_entries = species_filter(
        database_entries,
        mol_entries_pickle_location=folder + "/mol_entries.pickle",
        species_report=folder + "/unfiltered_species_report.tex",
        species_decision_tree=species_decision_tree,
        coordimer_weight=lambda mol: (mol.get_free_energy(params["temperature"])),
        species_logging_decision_tree=species_decision_tree,
        generate_unfiltered_mol_pictures=False,
    )

    print(len(mol_entries), "initial mol entries")

    bucket(mol_entries, folder + "/buckets.sqlite")

    print(len(mol_entries), "final mol entries")

    dispatcher_payload = DispatcherPayload(
        folder + "/buckets.sqlite",
        folder + "/rn.sqlite",
        folder + "/reaction_report.tex",
    )

    worker_payload = WorkerPayload(
        folder + "/buckets.sqlite",
        euvl_phase2_reaction_decision_tree,
        params,
        euvl_phase2_reaction_decision_tree,
    )

    dumpfn(dispatcher_payload, folder + "/dispatcher_payload.json")
    dumpfn(worker_payload, folder + "/worker_payload.json")

    subprocess.run(
        [
            "mpirun",
            "--use-hwthread-cpus",
            "-n",
            number_of_threads,
            "python",
            "run_network_generation.py",
            folder + "/mol_entries.pickle",
            folder + "/dispatcher_payload.json",
            folder + "/worker_payload.json",
        ]
    )

    phase1_network_loader = NetworkLoader(
        phase1_folder + "/rn.sqlite",
        phase1_folder + "/mol_entries.pickle",
        phase1_folder + f"/initial_state.sqlite",
    )
    phase1_network_loader.load_initial_state_and_trajectories()
    phase1_simulation_replayer = SimulationReplayer(phase1_network_loader)
    phase1_simulation_replayer.compute_trajectory_final_states()

    for seed in range(1000, 2000):

        initial_state = {}
        for ii, val in enumerate(phase1_simulation_replayer.final_states[seed]):
            if int(val) > 0:
                initial_state[ii] = int(val)

        insert_initial_state(
            initial_state,
            mol_entries,
            folder + "/initial_state_" + str(seed) + ".sqlite",
        )

        subprocess.run(
            [
                "GMC",
                "--reaction_database=" + folder + "/rn.sqlite",
                "--initial_state_database=" + folder + "/initial_state_" + str(seed) + ".sqlite",
                "--number_of_simulations=" + number_of_threads,
                "--base_seed=" + str(1000+(seed-1000)*int(number_of_threads)),
                "--thread_count=" + number_of_threads,
                "--step_cutoff=500",
            ]
        )

    network_loader = NetworkLoader(
        folder + "/rn.sqlite",
    folder + "/mol_entries.pickle",
    )

    for seed in range(1000, 2000):
        network_loader.set_initial_state_db(folder + "/initial_state_"+str(seed)+".sqlite")
        network_loader.load_initial_state_and_trajectories()

    report_generator = ReportGenerator(
        network_loader.mol_entries, folder + "/dummy.tex", rebuild_mol_pictures=True
    )
    reaction_tally_report(network_loader, folder + "/reaction_tally.tex", cutoff=10)
    species_report(network_loader, folder + "/species_report.tex")
    simulation_replayer = SimulationReplayer(network_loader)
    final_state_report(simulation_replayer, folder + "/final_state_report.tex")

    sink_report(simulation_replayer, folder + "/sink_report.tex")

    tps_plus1_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/tps.xyz", 1)
    phs_0_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/phs.xyz", 0)
    tba_0_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/tba.xyz", 0)
    nf_minus1_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/nf.xyz", -1)

    phase2_important_species = [tps_plus1_id, phs_0_id, tba_0_id, nf_minus1_id]

    colors = list(mcolors.TABLEAU_COLORS.values())
    phase2_colorstyle_list = []
    for ii, species in enumerate(phase2_important_species):        
        phase2_colorstyle_list.append([colors[ii], "solid"])

    ii = 0
    for mol_id in simulation_replayer.sinks:
        if mol_id not in phase2_important_species:
            phase2_important_species.append(mol_id)
            phase2_colorstyle_list.append([colors[ii%len(colors)], "dashed"])
            ii += 1

    phase1_important_species = copy.deepcopy(phase2_important_species)
    phase1_important_species.append(len(mol_entries))

    phase1_colorstyle_list = copy.deepcopy(phase2_colorstyle_list)
    phase1_colorstyle_list.append(["black", "dotted"])

    phase1_simulation_replayer.time_series_graph(
        seeds=[i for i in range(1000,2000)],
        species_of_interest=phase1_important_species,
        path=os.path.join(folder,"phase1_time_series"),
        custom_y_max=36,
        custom_colorstyle_list=phase1_colorstyle_list
    )

    simulation_replayer.time_series_graph(
        seeds=[i for i in range(1000,1000+1000*int(number_of_threads))],
        species_of_interest=phase2_important_species,
        path=os.path.join(folder,"phase2_time_series"),
        custom_y_max=36,
        custom_colorstyle_list=phase2_colorstyle_list
    )

    tests_passed = True
    print("Number of species:", network_loader.number_of_species)
    if network_loader.number_of_species == 103:
        print(bcolors.PASS + "euvl_phase_2_test: correct number of species" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "euvl_phase_2_test: correct number of species" + bcolors.ENDC)
        tests_passed = False

    print("Number of reactions:", network_loader.number_of_reactions)
    if network_loader.number_of_reactions == 3912:
        print(bcolors.PASS + "euvl_phase_2_test: correct number of reactions" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "euvl_phase_2_test: correct number of reactions" + bcolors.ENDC)
        tests_passed = False

    return tests_passed


def euvl_bondnet_test():

    phase1_folder = "./euvl_phase1_test"
    folder = "./scratch/euvl_phase2_test"
    subprocess.run(["mkdir", folder])

    mol_json = "./data/euvl_test_set.json"
    database_entries = loadfn(mol_json)

    ## HY
    bondnet_test_json = "./scratch/euvl_phase2_test/reaction_networks_graphs"
    subprocess.run(["mkdir", bondnet_test_json])
    ##

    species_decision_tree = euvl_species_decision_tree

    params = {
        "temperature": ROOM_TEMP+200.0,
        "electron_free_energy": 0.0,
    }

    mol_entries, dgl_molecules_dict = species_filter(
        database_entries,
        mol_entries_pickle_location=folder + "/mol_entries.pickle",
        species_report=folder + "/unfiltered_species_report.tex",
        species_decision_tree=species_decision_tree,
        coordimer_weight=lambda mol: (mol.get_free_energy(params["temperature"])),
        species_logging_decision_tree=species_decision_tree,
        generate_unfiltered_mol_pictures=False,
    )

    print(len(mol_entries), "initial mol entries")

    bucket(mol_entries, folder + "/buckets.sqlite")

    print(len(mol_entries), "final mol entries")

    dispatcher_payload = DispatcherPayload(
        folder + "/buckets.sqlite",
        folder + "/rn.sqlite",
        folder + "/reaction_report.tex",
        bondnet_test_json + "/test.json"
    )

    worker_payload = WorkerPayload(
        folder + "/buckets.sqlite",
        euvl_phase2_reaction_decision_tree,
        params,
        euvl_phase2_reaction_decision_tree,
    )

    dumpfn(dispatcher_payload, folder + "/dispatcher_payload.json")
    dumpfn(worker_payload, folder + "/worker_payload.json")

    subprocess.run(
        [
            "mpirun",
            "--use-hwthread-cpus",
            "-n",
            number_of_threads,
            "python",
            "run_network_generation.py",
            folder + "/mol_entries.pickle",
            folder + "/dispatcher_payload.json",
            folder + "/worker_payload.json",
        ]
    )

    phase1_network_loader = NetworkLoader(
        phase1_folder + "/rn.sqlite",
        phase1_folder + "/mol_entries.pickle",
        phase1_folder + f"/initial_state.sqlite",
    )
    phase1_network_loader.load_initial_state_and_trajectories()
    phase1_simulation_replayer = SimulationReplayer(phase1_network_loader)
    phase1_simulation_replayer.compute_trajectory_final_states()

    for seed in range(1000, 2000):

        initial_state = {}
        for ii, val in enumerate(phase1_simulation_replayer.final_states[seed]):
            if int(val) > 0:
                initial_state[ii] = int(val)

        insert_initial_state(
            initial_state,
            mol_entries,
            folder + "/initial_state_" + str(seed) + ".sqlite",
        )

        subprocess.run(
            [
                "GMC",
                "--reaction_database=" + folder + "/rn.sqlite",
                "--initial_state_database=" + folder + "/initial_state_" + str(seed) + ".sqlite",
                "--number_of_simulations=" + number_of_threads,
                "--base_seed=" + str(1000+(seed-1000)*int(number_of_threads)),
                "--thread_count=" + number_of_threads,
                "--step_cutoff=500",
            ]
        )

    network_loader = NetworkLoader(
        folder + "/rn.sqlite",
    folder + "/mol_entries.pickle",
    )

    for seed in range(1000, 2000):
        network_loader.set_initial_state_db(folder + "/initial_state_"+str(seed)+".sqlite")
        network_loader.load_initial_state_and_trajectories()

    report_generator = ReportGenerator(
        network_loader.mol_entries, folder + "/dummy.tex", rebuild_mol_pictures=True
    )
    reaction_tally_report(network_loader, folder + "/reaction_tally.tex", cutoff=10)
    species_report(network_loader, folder + "/species_report.tex")
    simulation_replayer = SimulationReplayer(network_loader)
    final_state_report(simulation_replayer, folder + "/final_state_report.tex")

    sink_report(simulation_replayer, folder + "/sink_report.tex")

    tps_plus1_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/tps.xyz", 1)
    phs_0_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/phs.xyz", 0)
    tba_0_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/tba.xyz", 0)
    nf_minus1_id = find_mol_entry_from_xyz_and_charge(mol_entries, "./xyz_files/nf.xyz", -1)

    phase2_important_species = [tps_plus1_id, phs_0_id, tba_0_id, nf_minus1_id]

    colors = list(mcolors.TABLEAU_COLORS.values())
    phase2_colorstyle_list = []
    for ii, species in enumerate(phase2_important_species):        
        phase2_colorstyle_list.append([colors[ii], "solid"])

    ii = 0
    for mol_id in simulation_replayer.sinks:
        if mol_id not in phase2_important_species:
            phase2_important_species.append(mol_id)
            phase2_colorstyle_list.append([colors[ii%len(colors)], "dashed"])
            ii += 1

    phase1_important_species = copy.deepcopy(phase2_important_species)
    phase1_important_species.append(len(mol_entries))

    phase1_colorstyle_list = copy.deepcopy(phase2_colorstyle_list)
    phase1_colorstyle_list.append(["black", "dotted"])

    phase1_simulation_replayer.time_series_graph(
        seeds=[i for i in range(1000,2000)],
        species_of_interest=phase1_important_species,
        path=os.path.join(folder,"phase1_time_series"),
        custom_y_max=36,
        custom_colorstyle_list=phase1_colorstyle_list
    )

    simulation_replayer.time_series_graph(
        seeds=[i for i in range(1000,1000+1000*int(number_of_threads))],
        species_of_interest=phase2_important_species,
        path=os.path.join(folder,"phase2_time_series"),
        custom_y_max=36,
        custom_colorstyle_list=phase2_colorstyle_list
    )

    tests_passed = True
    print("Number of species:", network_loader.number_of_species)
    if network_loader.number_of_species == 103:
        print(bcolors.PASS + "euvl_phase_2_test: correct number of species" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "euvl_phase_2_test: correct number of species" + bcolors.ENDC)
        tests_passed = False

    print("Number of reactions:", network_loader.number_of_reactions)
    if network_loader.number_of_reactions == 3912:
        print(bcolors.PASS + "euvl_phase_2_test: correct number of reactions" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "euvl_phase_2_test: correct number of reactions" + bcolors.ENDC)
        tests_passed = False

    return tests_passed

tests = [
    # mg_test,
    # li_test,
    # flicho_test,
    # co2_test,
    # euvl_phase1_test,
    # euvl_phase2_test,
    euvl_bondnet_test
]

for test in tests:
    if not test():
        exit(1)
