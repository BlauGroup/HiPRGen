import sys
from monty.serialization import loadfn
from species_filter import *
from bucketing import *

if len(sys.argv) != 4:
    print("usage: python run_bucketing.py database_json_file mol_entries_pickle_location bucket_db_location ")
    quit()

database_json_file = sys.argv[1]
mol_entries_pickle_file = sys.argv[2]
bucket_db_file = sys.argv[3]

database_entries = loadfn(database_json_file)
mol_entries = species_filter(database_entries, mol_entries_pickle_file)
Bucket(mol_entries, bucket_db_file)
