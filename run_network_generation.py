import sys
from HiPRGen.reaction_filter import *
import pickle
from mpi4py import MPI
from HiPRGen.constants import *
from HiPRGen.reaction_questions import decision_tree_dict


# python run_network_generation.py mol_entries_pickle_file bucket_db_file rn_db_location generation_report_location decision_tree_name electron_free_energy
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

mol_entries_pickle_file = sys.argv[1]
bucket_db_file = sys.argv[2]
rn_db_file = sys.argv[3]
report_file = sys.argv[4]
decision_tree = sys.argv[5]
electron_free_energy = sys.argv[6]

with open(mol_entries_pickle_file, 'rb') as f:
    mol_entries = pickle.load(f)

if rank == DISPATCHER_RANK:
    dispatcher(mol_entries,
               bucket_db_file,
               rn_db_file,
               report_file)

else:
    worker(mol_entries,
           bucket_db_file,
           reaction_decision_tree=decision_tree_dict[decision_tree],
           params=params_dict[electron_free_energy]
           )
