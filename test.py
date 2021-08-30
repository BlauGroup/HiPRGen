import os
import sys
import subprocess
import sqlite3
import pickle
from HiPRGen.mc_analysis import *
from HiPRGen.network_loader import *
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge
from monty.serialization import loadfn
from HiPRGen.species_filter import *
from HiPRGen.bucketing import *
from HiPRGen.report_generator import *
from HiPRGen.initial_state import *



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

    # there is one non local part of species filtering: we consider two molecules
    # to be equivalent of they have the same covalent bonds and we choose one such
    # molecule in each coordimer class.
    species_decision_tree = li_ec_species_decision_tree
    mol_entries = species_filter(
        database_entries,
        folder + '/mol_entries.pickle',
        folder + '/unfiltered_species_report.tex',
        species_decision_tree,
        lambda mol: mol.solvation_free_energy
    )

    bucket(mol_entries, folder + '/buckets.sqlite')



    reaction_decision_tree = 'li_ec_reaction_decision_tree'
    logging_decision_tree = 'li_ec_redox_logging_decision_tree'
    electron_free_energy = "-1.4"
    subprocess.run([
        'mpiexec',
        '--use-hwthread-cpus',
        '-n',
        number_of_threads,
        'python',
        'run_network_generation.py',
        folder + '/mol_entries.pickle',
        folder + '/buckets.sqlite',
        folder + '/rn.sqlite',
        folder + '/reaction_report.tex',
        reaction_decision_tree,
        logging_decision_tree,
        electron_free_energy
    ])


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

    initial_state = {
        Li_plus_id : 30,
        EC_id : 30
    }

    insert_initial_state(initial_state, mol_entries, folder + '/rn.sqlite')

    subprocess.run([
        'RNMC',
        '--database=' + folder + '/rn.sqlite',
        '--number_of_simulations=1000',
        '--base_seed=1000',
        '--thread_count=' + number_of_threads,
        '--step_cutoff=200'])


    network_loader = NetworkLoader(
        folder + '/rn.sqlite',
        folder + '/mol_entries.pickle'
        )


    report_generator = ReportGenerator(
        network_loader.mol_entries,
        folder + '/dummy.tex',
        rebuild_mol_pictures=True)

    reaction_tally_report(
        network_loader,
        folder + '/reaction_tally.tex'
    )

    pathfinding = Pathfinding(network_loader)
    pathfinding.generate_pathway_report(
        LEDC_id,
        folder + '/LEDC_pathways.tex',
        sort_by_frequency=False
    )

    species_report(network_loader, folder + '/species_report.tex')
    sink_report(network_loader, folder + '/sink_report.tex')
    consumption_report(network_loader,
                       LEDC_id,
                       folder + '/LEDC_consumption_report.tex')


    tests_passed = True
    if network_loader.number_of_species == 197:
        print(bcolors.PASS +
              "li_test: correct number of species" +
              bcolors.ENDC)
    else:
        print(bcolors.FAIL +
              "li_test: correct number of species" +
              bcolors.ENDC)
        tests_passed = False



    if network_loader.number_of_reactions == 4989:
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
    species_decision_tree = mg_g2_species_decision_tree
    reaction_decision_tree = 'mg_g2_reaction_decision_tree'
    logging_decision_tree = 'standard_logging_decision_tree'
    electron_free_energy = "-2.06"

    database_entries = loadfn(mol_json)
    mol_entries = species_filter(
        database_entries,
        folder + '/mol_entries.pickle',
        folder + '/unfiltered_species_report.tex',
        species_decision_tree,
        lambda mol: mol.solvation_free_energy
    )

    bucket(mol_entries, folder + '/buckets.sqlite')

    subprocess.run([
        'mpiexec',
        '--use-hwthread-cpus',
        '-n',
        number_of_threads,
        'python',
        'run_network_generation.py',
        folder + '/mol_entries.pickle',
        folder + '/buckets.sqlite',
        folder + '/rn.sqlite',
        folder + '/reaction_report.tex',
        reaction_decision_tree,
        logging_decision_tree,
        electron_free_energy
    ])

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
        mg_g2_plus_plus_id : 30
    }

    insert_initial_state(initial_state, mol_entries, folder + '/rn.sqlite')

    subprocess.run([
        'RNMC',
        '--database=' + folder + '/rn.sqlite',
        '--number_of_simulations=1000',
        '--base_seed=1000',
        '--thread_count=' + number_of_threads,
        '--step_cutoff=200'])


    network_loader = NetworkLoader(
        folder + '/rn.sqlite',
        folder + '/mol_entries.pickle'
        )

    report_generator = ReportGenerator(
        network_loader.mol_entries,
        folder + '/dummy.tex',
        rebuild_mol_pictures=True)

    reaction_tally_report(
        network_loader,
        folder + '/reaction_tally.tex'
    )

    pathfinding = Pathfinding(network_loader)
    pathfinding.generate_pathway_report(
        c2h6_id,
        folder + '/C2H6_pathways.tex',
        sort_by_frequency=False
    )

    pathfinding.generate_pathway_report(
        c2h4_id,
        folder + '/C2H4_pathways.tex',
        sort_by_frequency=False
    )



    species_report(network_loader, folder + '/species_report.tex')

    tests_passed = True
    if network_loader.number_of_species == 90:
        print(bcolors.PASS +
              "mg_test: correct number of species" +
              bcolors.ENDC)
    else:
        print(bcolors.FAIL +
              "mg_test: correct number of species" +
              bcolors.ENDC)
        tests_passed = False



    if network_loader.number_of_reactions == 731:
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
        break
