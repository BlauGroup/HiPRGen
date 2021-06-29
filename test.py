import os
import sqlite3
import pickle
from mc_analysis import *
from network_loader import *

class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'

network_loader = NetworkLoader(
    './scratch/rn.sqlite',
    './scratch/ronald_mol_entries.pickle'
    )


def test_species_filter(network_loader):
    with open('./data/ronald_mol_entries.pickle', 'rb') as f:
        mol_entries_data = pickle.load(f)



    if network_loader.mol_entries == mol_entries_data:
        print(bcolors.PASS + "passed: test_species_filter" + bcolors.ENDC)
        return
    else:
        print(bcolors.FAIL + "failed: test_species_filter" + bcolors.ENDC)
        quit()

def test_reaction_gen(network_loader):

    if network_loader.number_of_reactions == 112060:
        print(bcolors.PASS + "passed: test_reaction_gen" + bcolors.ENDC)
        return
    else:
        print(bcolors.FAIL + "failed: test_reaction_gen" + bcolors.ENDC)
        quit()


def test_mc_pathfinding(network_loader):

    LEDC_id = find_mol_entry_from_xyz_and_charge(
        network_loader.mol_entries,
        './xyz_files/LEDC.xyz',
        0)

    pathfinding = Pathfinding(network_loader)
    pathfinding.generate_pathway_report(121, './scratch/LEDC_pathways.tex')



test_species_filter(network_loader)
test_reaction_gen(network_loader)
test_mc_pathfinding(network_loader)
