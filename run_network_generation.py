import sys
from HiPRGen.reaction_filter import *
from HiPRGen.reaction_questions import *
import pickle
from mpi4py import MPI
from HiPRGen.constants import *


# python run_network_generation.py mol_entries_pickle_file bucket_db_file rn_db_location generation_report_location
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

mol_entries_pickle_file = sys.argv[1]
bucket_db_file = sys.argv[2]
rn_db_file = sys.argv[3]
report_file = sys.argv[4]

with open(mol_entries_pickle_file, 'rb') as f:
    mol_entries = pickle.load(f)

standard_reaction_decision_tree = [
    (partial(dG_above_threshold, 0.5), Terminal.DISCARD),
    # redox branch
    (is_redox_reaction, [
        (too_many_reactants_or_products, Terminal.DISCARD),
        (dcharge_too_large, Terminal.DISCARD),
        (reactant_and_product_not_isomorphic, Terminal.DISCARD),
        (default_true, Terminal.KEEP)
    ]),
    (partial(star_count_diff_above_threshold, 4), Terminal.DISCARD),
    (reaction_is_covalent_decomposable, Terminal.DISCARD),
    (concerted_metal_coordination, Terminal.DISCARD),
    (concerted_metal_coordination_one_product, Terminal.DISCARD),
    (concerted_metal_coordination_one_reactant, Terminal.DISCARD),
    (metal_coordination_passthrough, Terminal.KEEP),
    (fragment_matching_found, Terminal.KEEP),
    (default_true, Terminal.DISCARD)
    ]

if rank == DISPATCHER_RANK:
    dispatcher(mol_entries,
               bucket_db_file,
               rn_db_file,
               report_file)

else:
    worker(mol_entries,
           bucket_db_file,
           params={
               'temperature' : ROOM_TEMP,
               'electron_free_energy' : -1.4
           }
           )
