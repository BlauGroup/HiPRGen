from mpi4py import MPI
from itertools import combinations
from report_generator import ReportGenerator
import sqlite3
from time import localtime, strftime
from reaction_questions import standard_reaction_decision_tree, standard_logging_decision_tree, run_decision_tree
from constants import *
from enum import Enum

"""
Phases 3 & 4 run in paralell using MPI

Phase 3: reaction gen and filtering
input: a bucket labeled by atom count
output: a list of reactions from that bucket
description: Loop through all possible reactions in the bucket and apply the decision tree. This will run in parallel over each bucket.

Phase 4: collating and indexing
input: all the outputs of phase 3 as they are generated
output: reaction network database
description: the worker processes from phase 3 are sending their reactions to this phase and it is writing them to DB as it gets them. We can ensure that duplicates don't get generated in phase 3 which means we don't need extra index tables on the db.

the code in this file is designed to run on a compute cluster using MPI.
"""


create_metadata_table = """
    CREATE TABLE metadata (
            number_of_species   INTEGER NOT NULL,
            number_of_reactions INTEGER NOT NULL,
            factor_zero         REAL NOT NULL,
            factor_two          REAL NOT NULL,
            factor_duplicate    REAL NOT NULL
    );
"""

insert_metadata = """
    INSERT INTO metadata VALUES (?, ?, ?, ?, ?)
"""

create_reactions_table = """
    CREATE TABLE reactions (
            reaction_id         INTEGER NOT NULL PRIMARY KEY,
            number_of_reactants INTEGER NOT NULL,
            number_of_products  INTEGER NOT NULL,
            reactant_1          INTEGER NOT NULL,
            reactant_2          INTEGER NOT NULL,
            product_1           INTEGER NOT NULL,
            product_2           INTEGER NOT NULL,
            rate                REAL NOT NULL,
            dG                  REAL NOT NULL
    );
"""

insert_reaction = """
    INSERT INTO reactions VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


# TODO: structure these global variables better
DISPATCHER_RANK = 0

# message tags

# sent by workers to the dispatcher once they have finished initializing
# only sent once
INITIALIZATION_FINISHED = 0

# sent by workers to the dispatcher to request a new table
SEND_ME_A_TABLE = 1

# sent by dispatcher to workers when delivering a new table
HERE_IS_A_TABLE = 2

# sent by workers to the dispatcher when reaction passes db decision tree
NEW_REACTION_DB = 3

# sent by workers to the dispatcher when reaction passes logging decision tree
NEW_REACTION_LOGGING = 4

# requests which are handled by the dispatcher
requests_to_handle = [
    SEND_ME_A_TABLE,
    NEW_REACTION_DB,
    NEW_REACTION_LOGGING
    ]

class WorkerState(Enum):
    INITIALIZING = 0
    RUNNING = 1
    FINISHED = 2


def log_message(string, verbose):
    if verbose:
        print(
            '[' + strftime('%H:%M:%S', localtime()) + ']',
            string)


def dispatcher(
        mol_entries,
        bucket_db,
        rn_db,
        generation_report_path,
        commit_freq=1000,
        factor_zero=1.0,
        factor_two=1.0,
        factor_duplicate=1.0,
        verbose=True
):

    comm = MPI.COMM_WORLD
    table_list = []
    bucket_con = sqlite3.connect(bucket_db)
    bucket_cur = bucket_con.cursor()

    res = bucket_cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    for name in res:
        table = name[0]
        table_list.append(table)

    log_message("creating reaction network db", verbose)
    rn_con = sqlite3.connect(rn_db)
    rn_cur = rn_con.cursor()
    rn_cur.execute(create_metadata_table)
    rn_cur.execute(create_reactions_table)
    rn_con.commit()

    log_message("initializing report generator", verbose)
    report_generator = ReportGenerator(
        mol_entries,
        generation_report_path)

    worker_states = {}

    worker_ranks = [i for i in range(comm.Get_size()) if i != DISPATCHER_RANK]

    for i in worker_ranks:
        worker_states[i] = WorkerState.INITIALIZING

    for i in worker_states:
        # block, waiting for workers to initialize
        comm.recv(source=i, tag=INITIALIZATION_FINISHED)
        worker_states[i] = WorkerState.RUNNING

    log_message("all workers running", verbose)

    reaction_index = 0

    log_message("handling requests", verbose)


    while True:
        if WorkerState.RUNNING not in worker_states.values():
            break

        status = MPI.Status()
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()
        rank = status.Get_source()

        if tag == SEND_ME_A_TABLE:
            if len(table_list) == 0:
                comm.send(None, dest=rank, tag=HERE_IS_A_TABLE)
                worker_states[rank] = WorkerState.FINISHED
            else:
                next_table = table_list.pop()
                comm.send(next_table, dest=rank, tag=HERE_IS_A_TABLE)
                log_message("dispatched " + next_table, verbose)


        elif tag == NEW_REACTION_DB:
            reaction = data
            rn_cur.execute(
                insert_reaction,
                (reaction_index,
                 reaction['number_of_reactants'],
                 reaction['number_of_products'],
                 reaction['reactants'][0],
                 reaction['reactants'][1],
                 reaction['products'][0],
                 reaction['products'][1],
                 reaction['rate'],
                 reaction['dG']
                 ))

            reaction_index += 1
            if reaction_index % commit_freq == 0:
                rn_con.commit()


        elif tag == NEW_REACTION_LOGGING:

            reaction = data[0]
            decision_path = data[1]

            report_generator.emit_verbatim(
                '\n'.join([str(f) for f in decision_path]))
            report_generator.emit_reaction(reaction)
            report_generator.emit_newline()



    log_message("finalzing database and generation report", verbose)
    rn_cur.execute(
        insert_metadata,
        (len(mol_entries) + 1,
         reaction_index + 1,
        factor_zero,
        factor_two,
        factor_duplicate)
    )


    report_generator.finished()
    rn_con.commit()
    bucket_con.close()
    rn_con.close()


def worker(
        mol_entries,
        bucket_db,
        reaction_decision_tree=standard_reaction_decision_tree,
        logging_decision_tree=standard_logging_decision_tree,
        params={
            'temperature' : ROOM_TEMP,
            'electron_free_energy' : -1.4
            }
):

    comm = MPI.COMM_WORLD
    con = sqlite3.connect(bucket_db)
    cur = con.cursor()


    comm.send(None, dest=DISPATCHER_RANK, tag=INITIALIZATION_FINISHED)

    while True:
        comm.send(None, dest=DISPATCHER_RANK, tag=SEND_ME_A_TABLE)
        table = comm.recv(source=DISPATCHER_RANK, tag=HERE_IS_A_TABLE)

        if table is None:
            break


        bucket = []
        res = cur.execute("SELECT * FROM " + table)
        for pair in res:
            bucket.append(pair)

        for (reactants, products) in combinations(bucket, 2):
            reaction = {
                'reactants' : reactants,
                'products' : products,
                'number_of_reactants' : len([i for i in reactants if i != -1]),
                'number_of_products' : len([i for i in products if i != -1])}


            reverse_reaction = {
                'reactants' : reaction['products'],
                'products' : reaction['reactants'],
                'number_of_reactants' : reaction['number_of_products'],
                'number_of_products' : reaction['number_of_reactants'],
                'reverse' : reaction
            }

            # reaction atom mapping is one of the most expensive operations we do
            # it takes ~0.02 seconds. If we compute the atom mapping for a reaction
            # we don't need to also compute if for the reverse reaction, so we couple
            # reaction / reverse pairs to facilitate that.

            # this attribute is only here for performance reasons. Question functions
            # should absolutely not be touching it unless they are about to compute
            # an atom mapping
            reaction['reverse'] = reverse_reaction

            decision_pathway_forward = []
            decision_pathway_reverse = []
            if run_decision_tree(reaction,
                                 mol_entries,
                                 params,
                                 reaction_decision_tree,
                                 decision_pathway_forward
                                 ):

                comm.send(
                    reaction,
                    dest=DISPATCHER_RANK,
                    tag=NEW_REACTION_DB)


            if run_decision_tree(reverse_reaction,
                                 mol_entries,
                                 params,
                                 reaction_decision_tree,
                                 decision_pathway_reverse
                                 ):

                comm.send(
                    reverse_reaction,
                    dest=DISPATCHER_RANK,
                    tag=NEW_REACTION_DB)

            if run_decision_tree(reaction,
                                 mol_entries,
                                 params,
                                 logging_decision_tree):

                comm.send(
                    (reaction, decision_pathway_forward),
                    dest=DISPATCHER_RANK,
                    tag=NEW_REACTION_LOGGING)




            if run_decision_tree(reverse_reaction,
                                 mol_entries,
                                 params,
                                 logging_decision_tree):

                comm.send(
                    (reverse_reaction, decision_pathway_reverse),
                    dest=DISPATCHER_RANK,
                    tag=NEW_REACTION_LOGGING)
