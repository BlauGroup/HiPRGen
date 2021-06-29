import sys
from HiPRGen.initial_state import *
import pickle

if len(sys.argv) != 3:
    print("usage: python run_bucketing.py rn_database_file mol_entries_pickle_file")
    quit()

rn_db_file = sys.argv[1]
mol_entries_pickle = sys.argv[2]

with open(mol_entries_pickle, 'rb') as f:
    mol_entries = pickle.load(f)

Li_plus_id = find_mol_entry_from_xyz_and_charge(
    mol_entries,
    './xyz_files/Li.xyz',
    1)

EC_id = find_mol_entry_from_xyz_and_charge(
    mol_entries,
    './xyz_files/EC.xyz',
    0)

initial_state = {
    Li_plus_id : 30,
    EC_id : 30
}

insert_initial_state(initial_state, mol_entries, rn_db_file)

