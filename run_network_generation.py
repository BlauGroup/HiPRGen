import sys
import pickle
from mpi4py import MPI
from monty.serialization import loadfn
import multiprocessing as mp


# can't use fork system calls with older versions of MPI.
# somewhere either inside the python pickle library or the python sqlite library
# multiprocessing is used, which by default will fork.
# to avoid possible memory corruption, we force python to use spawn in these instances
mp.set_start_method('spawn')


from HiPRGen.reaction_filter import (
    dispatcher,
    worker,
    DISPATCHER_RANK
)


# python run_network_generation.py mol_entries_pickle_file dispatcher_payload.json worker_payload.json


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

mol_entries_pickle_file = sys.argv[1]
dispatcher_payload_json = sys.argv[2]
worker_payload_json = sys.argv[3]

with open(mol_entries_pickle_file, 'rb') as f:
    mol_entries = pickle.load(f)



if rank == DISPATCHER_RANK:
    dispatcher_payload = loadfn(dispatcher_payload_json)
    dispatcher(mol_entries,
               dispatcher_payload
               )

else:
    worker_payload = loadfn(worker_payload_json)
    worker(mol_entries,
           worker_payload
           )
