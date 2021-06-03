from monty.serialization import loadfn, dumpfn
from species_filter import *


def test_species_filter():
    ronald_LIBE = loadfn("./data/ronald_LIBE.json")
    ronald_mol_entries = species_filter(ronald_LIBE)
    ronald_mol_entries_from_disk = loadfn("./data/ronald_mol_entries.json")
    return ronald_mol_entries == ronald_mol_entries_from_disk


if test_species_filter():
    print("passed: test_species_filter")
else:
    print("failed: test_species_filter")
    quit()
