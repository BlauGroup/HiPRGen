import sys
from monty.serialization import loadfn, dumpfn
from species_filter import *
from bucketing import *
from reaction_gen import *

if len(sys.argv) != 4:
    print("usage: python run.py json_file bucket_db_location rn_db_location")
    quit()


database_entries = loadfn(sys.argv[1])
mol_entries = species_filter(database_entries)
Bucket(mol_entries, sys.argv[2])
dispatcher(mol_entries, sys.argv[2], sys.argv[3])
