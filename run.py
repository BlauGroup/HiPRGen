import sys
from monty.serialization import loadfn, dumpfn
from species_filter import *
from bucketing import *
from reaction_gen import *

if len(sys.argv) != 6:
    print("usage: python run.py database_json_file mol_entries_json bucket_db_location rn_db_location generation_report_location")
    quit()


database_entries = loadfn(sys.argv[1])
mol_entries = species_filter(database_entries, sys.argv[2])
Bucket(mol_entries, sys.argv[3])
dispatcher(mol_entries, sys.argv[3], sys.argv[4], sys.argv[5])
