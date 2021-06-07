from monty.serialization import loadfn, dumpfn
from species_filter import *
from bucketing import *
from reaction_gen import *

class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'



def test_species_filter(ronald_LIBE):
    ronald_mol_entries = species_filter(ronald_LIBE, None)
    ronald_mol_entries_from_disk = loadfn("./data/ronald_mol_entries.json")


    if ronald_mol_entries == ronald_mol_entries_from_disk:
        print(bcolors.PASS + "passed: test_species_filter" + bcolors.ENDC)
        return ronald_mol_entries
    else:
        print(bcolors.FAIL + "failed: test_species_filter" + bcolors.ENDC)
        quit()

def test_bucketing(mol_entries):
    b = Bucket(mol_entries, './scratch/buckets.sqlite')
    # TODO: write some test logic

def test_reaction_gen(mol_entries):
    d = dispatcher(mol_entries, './scratch/buckets.sqlite', './scratch/rn.sqlite', standard_decision_tree)
    # TODO: write some test logic

ronald_LIBE = loadfn("./data/ronald_LIBE.json")

ronald_mol_entries = test_species_filter(ronald_LIBE)
test_bucketing(ronald_mol_entries)
test_reaction_gen(ronald_mol_entries)
