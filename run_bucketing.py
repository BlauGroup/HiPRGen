import sys
from monty.serialization import loadfn
from species_filter import *
from bucketing import *
from report_generator import ReportGenerator


if len(sys.argv) != 5:
    print("usage: python run_bucketing.py database_json_file mol_entries_pickle_file bucket_db_file report_file")
    quit()

database_json_file = sys.argv[1]
mol_entries_pickle_file = sys.argv[2]
bucket_db_file = sys.argv[3]
report_file = sys.argv[4]

database_entries = loadfn(database_json_file)
mol_entries = species_filter(database_entries, mol_entries_pickle_file)
Bucket(mol_entries, bucket_db_file)

# generate a dummy report to get the mol pictures
report_generator = ReportGenerator(
    mol_entries,
    report_file
)
