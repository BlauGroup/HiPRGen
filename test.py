import os
import sys
import subprocess

from HiPRGen.network_loader import NetworkLoader
from monty.serialization import loadfn, dumpfn
from HiPRGen.species_filter import species_filter
from HiPRGen.bucketing import bucket
from HiPRGen.constants import ROOM_TEMP, Terminal
from HiPRGen.reaction_filter_payloads import (
    DispatcherPayload,
    WorkerPayload
)

from HiPRGen.species_questions import (
    mg_species_decision_tree,
    li_species_decision_tree,
)

from HiPRGen.reaction_questions import (
    default_reaction_decision_tree,

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
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'

# HiPRGen is organized as a pipeline, where all the relevent data is
# stored in a sqlite database between phases. For this reason, during
# a single run of the full pipeline, it makes sense to store all the
# relevent files in a single directory. We have two test sets, a lithium
# set and a magnesium set. Since the lithium test set is older, we shall
# document that instead of the mg test set.


if os.path.isdir('./scratch'):
    subprocess.run(['rm', '-r', './scratch'])

subprocess.run(['mkdir', './scratch'])


def li_test():

    # folder is the where we store all our intermediate databases
    folder = './scratch/li_test'
    subprocess.run(['mkdir', folder])

    # The initial input to the pipeline is a list of LIBE or MADEIRA
    # dataset entries. We provide two examples in the data foloder.
    mol_json = './data/ronald_LIBE.json'
    database_entries = loadfn(mol_json)
    # The first step of the HiPRGen pipeline is passing the input molecules
    # through the species decision tree to discard molecules.
    species_decision_tree = li_species_decision_tree

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
        mol_entries_pickle_location=folder + '/mol_entries.pickle',
        species_report=folder + '/unfiltered_species_report.tex',
        species_decision_tree=species_decision_tree,
        coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy),
    )

    # Once we have generated our molecule list, we generate the bucket database
    # which is how we break up the reaction filtering amongst all avaliable
    # workers. It gets stored in the buckets.sqlite database.
    bucket(mol_entries, folder + '/buckets.sqlite')

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

    params = {
        'temperature': ROOM_TEMP,
        'electron_free_energy': -1.4
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

    # The dispatcher and worker payloads are passed through the MPI barrier
    # as JSON blobs dispatcher_payload and worker_payload
    dumpfn(dispatcher_payload, folder + '/dispatcher_payload.json')
    dumpfn(worker_payload, folder + '/worker_payload.json')

    subprocess.run(
        [
            'mpirun',
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

    # The network loader builds a python object around a reaction network
    # and the molecules to make it easier to use them.
    network_loader = NetworkLoader(
        folder + '/rn.sqlite',
        folder + '/mol_entries.pickle'
        )

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
    subprocess.run(['mkdir', folder])

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
            'temperature': ROOM_TEMP,
            'electron_free_energy': -2.06
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

    network_loader = NetworkLoader(
        folder + '/rn.sqlite',
        folder + '/mol_entries.pickle',
        )

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
    # flicho_test
]

for test in tests:
    if not test():
        exit(1)
