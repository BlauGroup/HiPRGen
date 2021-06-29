import sys
from reaction_filter import *
import pickle
from mpi4py import MPI


# python run_network_generation.py mol_entries_pickle_file bucket_db_file rn_db_location generation_report_location
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

mol_entries_pickle_file = sys.argv[1]
bucket_db_file = sys.argv[2]
rn_db_file = sys.argv[3]
report_file = sys.argv[4]

with open(mol_entries_pickle_file, 'rb') as f:
    mol_entries = pickle.load(f)

if rank == DISPATCHER_RANK:
    dispatcher(mol_entries,
               bucket_db_file,
               rn_db_file,
               report_file)

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


else:
    worker(mol_entries,
           bucket_db_file,
           logging_decision_tree=standard_reaction_decision_tree
           )


