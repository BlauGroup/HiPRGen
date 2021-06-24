import os
import sqlite3
import pickle

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
        return
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

    if metadata == [(38, 156, 1.0, 1.0, 1.0)]:
        print(bcolors.PASS + "passed: test_reaction_gen" + bcolors.ENDC)
        return
    else:
        print(bcolors.FAIL + "failed: test_reaction_gen" + bcolors.ENDC)
        quit()


test_species_filter()
test_bucketing()
test_reaction_gen()
