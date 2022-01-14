import os
import sys
import subprocess
import sqlite3
import pickle


from HiPRGen.network_loader import NetworkLoader
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge
from monty.serialization import loadfn, dumpfn
from HiPRGen.species_filter import species_filter
from HiPRGen.bucketing import bucket
from HiPRGen.report_generator import ReportGenerator
from HiPRGen.initial_state import insert_initial_state
from HiPRGen.constants import ROOM_TEMP, Terminal
from HiPRGen.reaction_filter_payloads import (
    DispatcherPayload,
    WorkerPayload
)

from HiPRGen.species_questions import (
    mg_species_decision_tree,
    li_species_decision_tree,
    positive_penalty,
    species_default_true
)

from HiPRGen.reaction_questions import (
    default_reaction_decision_tree,

)

from HiPRGen.mc_analysis import (
    reaction_tally_report,
    species_report,
    Pathfinding,
    SimulationReplayer,
    generate_pathway_report,
    sink_report,
    consumption_report,
    redox_report
)

# Since HiPRGen uses an end-to-end testing approach rather than testing
# each individual function, we have decided to use the tests as
# documentation, by explaining every single line.


# the first thing you need to think about when using HiPRGen is how many
# worker threads do you want to run. HiPRGen can be run with a single
# thread or thousands distrubuted across several nodes. For reaction
# networks with between 5000 and 10000 species, we have found that the
# optimal number of worker threads is between 1000 and 2000. If you try
# and use more than that, the worker threads are going to spend lots of
# time waiting for the dispatcher to get through all of the reactions it
# is being sent, which slows everything down. Fixing this would require
# a more complex distrubuted system, but it hasn't been an issue yet for
# the large reaction networks we have been generating.
if len(sys.argv) != 2:
    print("usage: python test.py number_of_threads")
    quit()


number_of_threads = sys.argv[1]

class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'


# HiPRGen is organized as a pipeline, where all the relevent data is
# stored in a sqlite database inbetween phases. For this reason, during
# a single run of the full pipeline, it makes sense to store all the
# relevent files in a single directory. We have two test sets, a lithium
# set and a magnesium set. Since the lithium test set is older, we shall
# document that instead of the mg test set.

if os.path.isdir('./scratch'):
    subprocess.run(['rm', '-r', './scratch'])

subprocess.run(['mkdir', './scratch'])
subprocess.run(['mkdir', './scratch/li_test' ])
subprocess.run(['mkdir', './scratch/mg_test' ])


def li_test():


    # folder is the where we store all our intermediate databases
    folder = './scratch/li_test'

    # the initial input to the pipeline is a list of LIBE or MADEIRA
    # dataset entries. We provide two examples in the data foloder.
    mol_json = './data/ronald_LIBE.json'
    database_entries = loadfn(mol_json)
    # the first step of the HiPRGen pipeline is passing the input molecules
    # through the a species decision tree to discard molecules. This happens
    # here rather than further complicating the DFT pipelines which generate the
    # input data for HiPRGen.
    species_decision_tree = li_species_decision_tree


    # there is one non local part of species filtering: we consider two
    # molecules to be equivalent of they have the same covalent bonds
    # and we choose one such molecule in each coordimer class using the
    # coodimer weight function.  Since most of our logging later on is
    # defined in terms of a fixed molecule set, logging for the species
    # filtering phase is messy, so ignore the species_report argument for
    # now. The second argument is where we store a pickle of the
    # filtered molecule entries for use in later phases.

    mol_entries = species_filter(
        database_entries,
        mol_entries_pickle_location=folder + '/mol_entries.pickle',
        species_report=folder + '/unfiltered_species_report.tex',
        species_decision_tree=species_decision_tree,
        coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy),
        species_logging_decision_tree=[
            (positive_penalty(), Terminal.KEEP),
            (species_default_true(), Terminal.DISCARD)],
        generate_unfiltered_mol_pictures=True
    )


    # once we have generated our molecule list, we generate the bucket database
    # which is how we break up the reaction filtering amoungst all avaliable workers.
    # it gets stored in the buckets.sqlite database.
    bucket(mol_entries, folder + '/buckets.sqlite')


    # reaction filtering is paralellized using MPI, so we need to spawn
    # an MPI instance to run it. That is why we can't just start
    # reaction filtering by calling a python function. We pass the
    # reaction decision tree, the logging decision tree and the electron
    # free energy as strings across this barrier. Every possible
    # reaction gets passed through both the reaction decision tree and
    # the logging decision tree.  if a reaction passes the reaction
    # decision tree, it gets written to the network.  if a reaction
    # passes the logging decision tree, it gets logged to the reaction
    # report along with what happened to it in reaction_decision_tree.

    # the reaction decision trees are constructed in
    # HiPRGen.reaction_questions

    params = {
        'temperature' : ROOM_TEMP,
        'electron_free_energy' : -1.4
    }

    dispatcher_payload = DispatcherPayload(
        folder + '/buckets.sqlite',
        folder + '/rn.sqlite',
        folder + '/reaction_report.tex'
    )

    worker_payload = WorkerPayload(
        folder + '/buckets.sqlite',
        default_reaction_decision_tree,
        params,
        Terminal.DISCARD
    )


    # the dispatcher and worker payloads are passed through the MPI barrier
    # as json blobs dispatcher_payload and worker_payload
    dumpfn(dispatcher_payload, folder + '/dispatcher_payload.json')
    dumpfn(worker_payload, folder + '/worker_payload.json')

    subprocess.run(
        [
            'mpiexec',
            '--use-hwthread-cpus',
            '-n',
            number_of_threads,
            'python',
            'run_network_generation.py',
            folder + '/mol_entries.pickle',
            folder + '/dispatcher_payload.json',
            folder + '/worker_payload.json'
        ]
    )

    # after we have generated the mol_entries, we refer to molecules by
    # their index. The function find_mol_entry_from_xyz_and_charge is able
    # to find a mol entry just from the xyz positions of its atoms, although
    # it isn't 100% reliable.
    Li_plus_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        './xyz_files/Li.xyz',
        1)

    EC_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        './xyz_files/EC.xyz',
        0)

    LEDC_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        './xyz_files/LEDC.xyz',
        0)


    # after generating a reaction network, it is stored in rn.sqlite. We want
    # to use Monte Carlo simulation to better understand the network, and for that
    # we need to insert an initial condition.
    initial_state = {
        Li_plus_id : 30,
        EC_id : 30
    }

    # the initial state and trajectories (after simulation) are stored in
    # a seperate database, in this case called initial_state.sqlite.
    # This facilitates running multiple simulations of the same network
    # with different initial conditions at the same time.
    insert_initial_state(initial_state, mol_entries, folder + '/initial_state.sqlite')


    # GMC is a high performance reaction network monte carlo simulator using the
    # gillespie algorithm: https://github.com/BlauGroup/RNMC
    subprocess.run([
        'GMC',
        '--reaction_database=' + folder + '/rn.sqlite',
        '--initial_state_database=' + folder + '/initial_state.sqlite',
        '--number_of_simulations=1000',
        '--base_seed=1000',
        '--thread_count=' + number_of_threads,
        '--step_cutoff=200',
        '--dependency_threshold=1'
    ])

    # the network loader builds a python object around a reaction network
    # and the molecules to make it easier to use them.
    network_loader = NetworkLoader(
        folder + '/rn.sqlite',
        folder + '/mol_entries.pickle',
        folder + '/initial_state.sqlite'
        )

    network_loader.load_trajectories()
    network_loader.load_initial_state()



    # HiPRGen has analysis tools to understand what happened in our simulation.
    # the output files are written into the same folder as where the reaction
    # network is stored.

    # this report is empty, but we use it to generate the molecule pictures.
    # this is an expensive operation, so we only want do do it once.
    report_generator = ReportGenerator(
        network_loader.mol_entries,
        folder + '/dummy.tex',
        rebuild_mol_pictures=True)

    reaction_tally_report(
        network_loader,
        folder + '/reaction_tally.tex'
    )

    # pathfinding is the main goal of HiPRGen / GMC.
    # run pdflatex LEDC_pathways.tex to see all the ways that LEDC was
    # produced in the simulations of our lithium test network. Note that this
    # network has ~5000 reactions. Our production networks have
    # between 50-100 million reactions.

    redox_report(network_loader, folder + '/redox_report.tex', params)

    pathfinding = Pathfinding(network_loader)
    generate_pathway_report(
        pathfinding,
        LEDC_id,
        folder + '/LEDC_pathways.tex',
        sort_by_frequency=False
    )

    generate_pathway_report(
        pathfinding,
        130,
        folder + '/130_pathways.tex',
        sort_by_frequency=False
    )



    species_report(network_loader, folder + '/species_report.tex')


    simulation_replayer = SimulationReplayer(network_loader)
    sink_report(simulation_replayer, folder + '/sink_report.tex')
    consumption_report(simulation_replayer,
                       LEDC_id,
                       folder + '/LEDC_consumption_report.tex')


    tests_passed = True
    if network_loader.number_of_species == 190:
        print(bcolors.PASS +
              "li_test: correct number of species" +
              bcolors.ENDC)
    else:
        print(bcolors.FAIL +
              "li_test: correct number of species" +
              bcolors.ENDC)
        tests_passed = False



    if network_loader.number_of_reactions == 4921:
        print(bcolors.PASS +
              "li_test: correct number of reactions" +
              bcolors.ENDC)
    else:
        print(bcolors.FAIL +
              "li_test: correct number of reactions" +
              bcolors.ENDC)
        tests_passed = False

    return tests_passed



def mg_test():


    folder = './scratch/mg_test'
    mol_json = './data/sam_G2.json'
    species_decision_tree = mg_species_decision_tree

    database_entries = loadfn(mol_json)



    mol_entries = species_filter(
        database_entries,
        folder + '/mol_entries.pickle',
        folder + '/unfiltered_species_report.tex',
        species_decision_tree,
        coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy)
    )



    bucket(mol_entries, folder + '/buckets.sqlite')


    dispatcher_payload = DispatcherPayload(
        folder + '/buckets.sqlite',
        folder + '/rn.sqlite',
        folder + '/reaction_report.tex'
    )

    worker_payload = WorkerPayload(
        folder + '/buckets.sqlite',
        default_reaction_decision_tree,
        {
            'temperature' : ROOM_TEMP,
            'electron_free_energy' : -2.06
        },
        Terminal.DISCARD
    )


    dumpfn(dispatcher_payload, folder + '/dispatcher_payload.json')
    dumpfn(worker_payload, folder + '/worker_payload.json')

    subprocess.run(
        [
            'mpiexec',
            '--use-hwthread-cpus',
            '-n',
            number_of_threads,
            'python',
            'run_network_generation.py',
            folder + '/mol_entries.pickle',
            folder + '/dispatcher_payload.json',
            folder + '/worker_payload.json'
        ]
    )


    mg_g2_plus_plus_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        './xyz_files/mgg2.xyz',
        2)

    c2h4_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        './xyz_files/c2h4.xyz',
        0)

    c2h6_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        './xyz_files/c2h6.xyz',
        0)

    initial_state = {
        33 : 30,
        81 : 30
    }


    insert_initial_state(initial_state, mol_entries, folder + '/initial_state.sqlite')


    subprocess.run([
        'GMC',
        '--reaction_database=' + folder + '/rn.sqlite',
        '--initial_state_database=' + folder + '/initial_state.sqlite',
        '--number_of_simulations=1000',
        '--base_seed=1000',
        '--thread_count=' + number_of_threads,
        '--step_cutoff=200',
        '--dependency_threshold=1'
    ])



    network_loader = NetworkLoader(
        folder + '/rn.sqlite',
        folder + '/mol_entries.pickle',
        folder + '/initial_state.sqlite'
        )

    network_loader.load_trajectories()
    network_loader.load_initial_state()



    report_generator = ReportGenerator(
        network_loader.mol_entries,
        folder + '/dummy.tex',
        rebuild_mol_pictures=True)

    reaction_tally_report(
        network_loader,
        folder + '/reaction_tally.tex'
    )

    pathfinding = Pathfinding(network_loader)

    generate_pathway_report(
        pathfinding,
        c2h6_id,
        folder + '/C2H6_pathways.tex',
        sort_by_frequency=False
    )

    generate_pathway_report(
        pathfinding,
        c2h4_id,
        folder + '/C2H4_pathways.tex',
        sort_by_frequency=False
    )



    species_report(network_loader, folder + '/species_report.tex')

    tests_passed = True
    if network_loader.number_of_species == 83:
        print(bcolors.PASS +
              "mg_test: correct number of species" +
              bcolors.ENDC)
    else:
        print(bcolors.FAIL +
              "mg_test: correct number of species" +
              bcolors.ENDC)
        tests_passed = False



    if network_loader.number_of_reactions == 788:
        print(bcolors.PASS +
              "mg_test: correct number of reactions" +
              bcolors.ENDC)
    else:
        print(bcolors.FAIL +
              "mg_test: correct number of reactions" +
              bcolors.ENDC)
        tests_passed = False

    return tests_passed


tests = [
    mg_test,
    li_test,
]

for test in tests:
    if not test():
        exit(1)

