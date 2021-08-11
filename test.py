import os
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


class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'

if os.path.isdir('./scratch'):
    subprocess.run(['rm', '-r', './scratch'])

subprocess.run(['mkdir', './scratch'])
subprocess.run(['mkdir', './scratch/li_test' ])



def li_test():
    folder = './scratch/li_test'
    database_entries = loadfn('./data/ronald_LIBE.json')
    mol_entries = species_filter(
        database_entries,
        folder + '/ronald_mol_entries.pickle',
        folder + '/report_file',
        li_species_decision_tree
    )

    bucket(mol_entries, folder + '/buckets.sqlite')

    subprocess.run([
        'mpiexec',
        '--use-hwthread-cpus',
        '-n',
        '8',
        'python',
        'run_network_generation.py',
        folder + '/ronald_mol_entries.pickle',
        folder + '/buckets.sqlite',
        folder + '/rn.sqlite',
        folder + '/reaction_report.tex',
        'standard_reaction_decision_tree'
    ])


    Li_plus_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        './xyz_files/Li.xyz',
        1)

    EC_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        './xyz_files/EC.xyz',
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
        '--thread_count=8',
        '--step_cutoff=200'])


    network_loader = NetworkLoader(
        folder + '/rn.sqlite',
        folder + '/ronald_mol_entries.pickle'
        )


    report_generator = ReportGenerator(
        network_loader.mol_entries,
        folder + '/dummy.tex',
        rebuild_mol_pictures=True)

    reaction_tally_report(
        network_loader,
        folder + '/reaction_tally.tex'
    )

    LEDC_id = find_mol_entry_from_xyz_and_charge(
        network_loader.mol_entries,
        './xyz_files/LEDC.xyz',
        0)

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


    if network_loader.number_of_species == 197:
        print(bcolors.PASS +
              "li_test: correct number of species" +
              bcolors.ENDC)
    else:
        print(bcolors.FAIL +
              "li_test: correct number of species" +
              bcolors.ENDC)



    if network_loader.number_of_reactions == 5525:
        print(bcolors.PASS +
              "li_test: correct number of reactions" +
              bcolors.ENDC)
    else:
        print(bcolors.FAIL +
              "li_test: correct number of reactions" +
              bcolors.ENDC)

li_test()
