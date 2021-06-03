from monty.serialization import loadfn, dumpfn
from species_filter import *

class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'



def test_species_filter(ronald_LIBE):
    ronald_mol_entries = species_filter(ronald_LIBE, None)
    ronald_mol_entries_from_disk = loadfn("./data/ronald_mol_entries.json")
    if ronald_mol_entries == ronald_mol_entries_from_disk:
        print(bcolors.PASS + "passed: test_species_filter" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "failed: test_species_filter" + bcolors.ENDC)
        quit()

    return ronald_mol_entries



ronald_LIBE = loadfn("./data/ronald_LIBE.json")

ronald_mol_entries = test_species_filter(ronald_LIBE)
