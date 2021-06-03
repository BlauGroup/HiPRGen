import sys
from monty.serialization import loadfn, dumpfn
from species_filter import *
from bucketing import *

if len(sys.argv) != 3:
    print("usage: python run.py json_file db_location")
    quit()


database_entries = loadfn(sys.argv[1])
mol_entries = species_filter(database_entries, None)
Bucket(mol_entries, sys.argv[2])
