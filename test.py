from monty.serialization import loadfn, dumpfn
from species_filter import *
from bucketing import *
from reaction_gen import *
from functools import partial
import os
import sqlite3

class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'

if os.path.isdir('./scratch'):
    os.system('rm -r ./scratch')

os.system('mkdir scratch')



def test_species_filter(ronald_LIBE):
    ronald_mol_entries = species_filter(
        ronald_LIBE,
        "./scratch/ronald_mol_entries.pickle",
        verbose=False)

    ronald_mol_entries_from_disk = loadfn("./data/ronald_mol_entries.json")


    if ronald_mol_entries == ronald_mol_entries_from_disk:
        print(bcolors.PASS + "passed: test_species_filter" + bcolors.ENDC)
        return ronald_mol_entries
    else:
        print(bcolors.FAIL + "failed: test_species_filter" + bcolors.ENDC)
        quit()

def test_bucketing(mol_entries):
    b = Bucket(mol_entries, './scratch/buckets.sqlite')
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




def test_reaction_gen(mol_entries):
    d = dispatcher(
        mol_entries,
        './scratch/buckets.sqlite',
        './scratch/rn.sqlite',
        './scratch/generation_report.tex',
        verbose=False
    )

    con = sqlite3.connect('./scratch/rn.sqlite')
    cur = con.cursor()
    metadata = list(cur.execute("SELECT * FROM metadata"))

    if metadata == [(38, 156, 1.0, 1.0, 1.0)]:
        print(bcolors.PASS + "passed: test_reaction_gen" + bcolors.ENDC)
        return
    else:
        print(bcolors.FAIL + "failed: test_reaction_gen" + bcolors.ENDC)
        quit()



ronald_LIBE = loadfn("./data/ronald_LIBE.json")

ronald_mol_entries = test_species_filter(ronald_LIBE)
test_bucketing(ronald_mol_entries)
test_reaction_gen(ronald_mol_entries)


