import os
import sqlite3
import pickle
from initial_state import *
from mc_analysis import *

class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'


def test_species_filter():
    with open('./scratch/ronald_mol_entries.pickle', 'rb') as f:
        mol_entries_scratch = pickle.load(f)

    with open('./data/ronald_mol_entries.pickle', 'rb') as f:
        mol_entries_data = pickle.load(f)



    if mol_entries_scratch == mol_entries_data:
        print(bcolors.PASS + "passed: test_species_filter" + bcolors.ENDC)
        return mol_entries_scratch
    else:
        print(bcolors.FAIL + "failed: test_species_filter" + bcolors.ENDC)
        quit()

def test_bucketing():
    con = sqlite3.connect('./scratch/buckets.sqlite')
    cur = con.cursor()
    number_of_buckets = list(cur.execute(
        "SELECT count(*) FROM sqlite_master WHERE type = 'table'"))[0][0]


    if number_of_buckets == 178:
        print(bcolors.PASS + "passed: test_bucketing" + bcolors.ENDC)
        return
    else:
        print(bcolors.FAIL + "failed: test_bucketing" + bcolors.ENDC)
        quit()

def test_reaction_gen():
    con = sqlite3.connect('./scratch/rn.sqlite')
    cur = con.cursor()
    metadata = list(cur.execute("SELECT * FROM metadata"))

    if metadata == [(302, 112060, 1.0, 1.0, 1.0)]:
        print(bcolors.PASS + "passed: test_reaction_gen" + bcolors.ENDC)
        return
    else:
        print(bcolors.FAIL + "failed: test_reaction_gen" + bcolors.ENDC)
        quit()


def test_initial_state(mol_entries):
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

    insert_initial_state(initial_state, mol_entries, './scratch/rn.sqlite')

    # there isn't really anything to test here. Just making sure the code runs
    print(bcolors.PASS + "passed: test_initial_state" + bcolors.ENDC)


def test_mc_pathfinding(mol_entries):

    LEDC_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        './xyz_files/LEDC.xyz',
        0)

    pathfinding = Pathfinding(
        './scratch/rn.sqlite',
        mol_entries,
        './scratch/LEDC_pathways.tex'
    )

    pathfinding.generate_pathway_report(LEDC_id)



mol_entries = test_species_filter()
test_bucketing()
test_reaction_gen()
# test_initial_state(mol_entries)
test_mc_pathfinding(mol_entries)

