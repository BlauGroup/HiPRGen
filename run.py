import sys
from monty.serialization import loadfn, dumpfn
from species_filter import *
from reaction_filter import *
from bucketing import *
import pickle

if len(sys.argv) != 6:
    print("usage: python run.py database_json_file mol_entries_pickle_location bucket_db_location rn_db_location generation_report_location")
    quit()


database_entries = loadfn(sys.argv[1])
mol_entries = species_filter(database_entries, sys.argv[2])

# with open(sys.argv[2], 'rb') as f:
#     mol_entries = pickle.load(f)

Bucket(mol_entries, sys.argv[3])
dispatcher(mol_entries, sys.argv[3], sys.argv[4], sys.argv[5])
