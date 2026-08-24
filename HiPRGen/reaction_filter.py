from mpi4py import MPI
from itertools import permutations, product
from HiPRGen.report_generator import ReportGenerator
from time import time
from HiPRGen.logging import log_message
import sqlite3
from enum import Enum
from math import floor

from HiPRGen.reaction_questions import (
    run_decision_tree
)

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
            number_of_reactions INTEGER NOT NULL
    );
"""

insert_metadata = """
    INSERT INTO metadata VALUES (?, ?)
"""

# it is important that reaction_id is the primary key
# otherwise the network loader will be extremely slow.
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
            dG                  REAL NOT NULL,
            dG_barrier          REAL NOT NULL,
            is_redox            INTEGER NOT NULL
    );
"""


insert_reaction = """
    INSERT INTO reactions VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""

get_complex_group_sql = """
    SELECT * FROM complexes WHERE composition_id=? AND group_id=?
"""


# TODO: structure these global variables better
DISPATCHER_RANK = 0

# message tags

# sent by workers to the dispatcher once they have finished initializing
# only sent once
INITIALIZATION_FINISHED = 0

# sent by workers to the dispatcher to request a new table
SEND_ME_A_WORK_BATCH = 1

# sent by dispatcher to workers when delivering a new table
HERE_IS_A_WORK_BATCH = 2

# sent by workers to the dispatcher when reaction passes db decision tree
NEW_REACTION_DB = 3

# sent by workers to the dispatcher when reaction passes logging decision tree
NEW_REACTION_LOGGING = 4

class WorkerState(Enum):
    INITIALIZING = 0
    RUNNING = 1
    FINISHED = 2

def dispatcher(
        mol_entries,
        dispatcher_payload
):

    comm = MPI.COMM_WORLD
    work_batch_list = []
    bucket_con = sqlite3.connect(dispatcher_payload.bucket_db_file)
    bucket_cur = bucket_con.cursor()
    size_cur = bucket_con.cursor()

    res = bucket_cur.execute("SELECT * FROM group_counts")
    for (composition_id, count) in res:
        for (i,j) in product(range(count), repeat=2):
            work_batch_list.append(
                (composition_id, i, j))

    # Sort ascending by estimated batch cost (the product of the two
    # groups' member counts) so that the main loop's
    # work_batch_list.pop() - which removes from the end of the list -
    # hands out the largest, slowest batches first. This is the
    # Longest-Processing-Time-first (LPT) scheduling heuristic:
    # front-loading the big batches and saving the smallest ones for last
    # keeps the run's tail short, since by the time only a few batches
    # remain they're all cheap, so no single straggler batch can leave
    # the rest of the worker fleet idle waiting on it.
    group_member_counts = {}
    for (composition_id, group_id, member_count) in bucket_cur.execute(
            "SELECT composition_id, group_id, COUNT(*) FROM complexes "
            "GROUP BY composition_id, group_id"):
        group_member_counts[(composition_id, group_id)] = member_count

    work_batch_list.sort(
        key=lambda batch: (
            group_member_counts[(batch[0], batch[1])] *
            group_member_counts[(batch[0], batch[2])]))

    composition_names = {}
    res = bucket_cur.execute("SELECT * FROM compositions")
    for (composition_id, composition) in res:
        composition_names[composition_id] = composition

    log_message("creating reaction network db")
    rn_con = sqlite3.connect(dispatcher_payload.reaction_network_db_file)
    # rn.sqlite is only ever opened here - a single connection from the
    # dispatcher, for the lifetime of the run - so there's no
    # multi-connection exposure for WAL's usual network-filesystem caveat
    # (its shared-memory wal-index doesn't coordinate reliably across
    # processes on Lustre/NFS) to apply to. WAL replaces the rollback
    # journal's per-commit create/delete-file churn with appends to a
    # single -wal file, and synchronous=NORMAL (safe under WAL - the WAL
    # file itself is the durability mechanism) skips an fsync per commit;
    # both matter when this is a high-commit-rate file on a
    # Lustre-backed path.
    rn_con.execute("PRAGMA journal_mode=WAL")
    rn_con.execute("PRAGMA synchronous=NORMAL")
    rn_cur = rn_con.cursor()
    rn_cur.execute(create_metadata_table)
    rn_cur.execute(create_reactions_table)
    rn_con.commit()

    log_message("initializing report generator")

    # since MPI processes spin lock, we don't want to have the dispathcer
    # spend a bunch of time generating molecule pictures
    report_generator = ReportGenerator(
        mol_entries,
        dispatcher_payload.report_file,
        rebuild_mol_pictures=False
    )

    worker_states = {}

    worker_ranks = [i for i in range(comm.Get_size()) if i != DISPATCHER_RANK]

    for i in worker_ranks:
        worker_states[i] = WorkerState.INITIALIZING

    for i in worker_states:
        # block, waiting for workers to initialize
        comm.recv(source=i, tag=INITIALIZATION_FINISHED)
        worker_states[i] = WorkerState.RUNNING

    log_message("all workers running")

    reaction_index = 0

    log_message("handling requests")

    batches_left_at_last_checkpoint = len(work_batch_list)
    last_checkpoint_time = floor(time())
    while True:
        if WorkerState.RUNNING not in worker_states.values():
            break

        current_time = floor(time())
        time_diff = current_time - last_checkpoint_time
        if ( current_time % dispatcher_payload.checkpoint_interval == 0 and
             time_diff > 0):
            batches_left_at_current_checkpoint = len(work_batch_list)
            batch_count_diff = (
                batches_left_at_last_checkpoint -
                batches_left_at_current_checkpoint)

            batch_consumption_rate = batch_count_diff / time_diff

            log_message("batches remaining:", batches_left_at_current_checkpoint)
            log_message("batch consumption rate:",
                        batch_consumption_rate,
                        "batches per second")


            batches_left_at_last_checkpoint = batches_left_at_current_checkpoint
            last_checkpoint_time = current_time


        status = MPI.Status()
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()
        rank = status.Get_source()

        if tag == SEND_ME_A_WORK_BATCH:
            if len(work_batch_list) == 0:
                comm.send(None, dest=rank, tag=HERE_IS_A_WORK_BATCH)
                worker_states[rank] = WorkerState.FINISHED
            else:
                # pop removes and returns the last item in the list
                work_batch = work_batch_list.pop()
                comm.send(work_batch, dest=rank, tag=HERE_IS_A_WORK_BATCH)
                composition_id, group_id_0, group_id_1 = work_batch
                log_message(
                    "dispatched",
                    composition_names[composition_id],
                    ": group ids:",
                    group_id_0, group_id_1
                )


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
                 reaction['dG'],
                 reaction['dG_barrier'],
                 reaction['is_redox']
                 ))

            reaction_index += 1
            if reaction_index % dispatcher_payload.commit_frequency == 0:
                rn_con.commit()


        elif tag == NEW_REACTION_LOGGING:

            reaction = data[0]
            decision_path = data[1]

            report_generator.emit_verbatim(decision_path)
            report_generator.emit_reaction(reaction)
            report_generator.emit_bond_breakage(reaction)
            report_generator.emit_newline()



    log_message("finalzing database and generation report")
    rn_cur.execute(
        insert_metadata,
        (len(mol_entries),
         reaction_index)
    )


    report_generator.finished()
    rn_con.commit()
    # Fold the WAL file back into rn.sqlite and drop back to the default
    # rollback-journal mode, so everything downstream (report generation,
    # NetworkLoader, ad hoc sqlite3 opens from analysis scripts) sees an
    # ordinary single-file db with no WAL/shm machinery to worry about.
    rn_con.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    rn_con.execute("PRAGMA journal_mode=DELETE")
    bucket_con.close()
    rn_con.close()


def worker(
        mol_entries,
        worker_payload
):

    comm = MPI.COMM_WORLD
    con = sqlite3.connect(worker_payload.bucket_db_file)
    cur = con.cursor()


    comm.send(None, dest=DISPATCHER_RANK, tag=INITIALIZATION_FINISHED)

    while True:
        comm.send(None, dest=DISPATCHER_RANK, tag=SEND_ME_A_WORK_BATCH)
        work_batch = comm.recv(source=DISPATCHER_RANK, tag=HERE_IS_A_WORK_BATCH)

        if work_batch is None:
            break

        # fragment_matching_found (HiPRGen.reaction_questions) memoizes
        # its per-complex work in this dict when present, since the same
        # complex recurs many times across a batch's pairings and its
        # per-complex result never depends on the pairing partner. Reset
        # fresh each batch so it never grows past what a single batch's
        # complexes need.
        worker_payload.params['_fragment_summary_cache'] = {}

        composition_id, group_id_0, group_id_1 = work_batch


        if group_id_0 == group_id_1:

            res = cur.execute(
                get_complex_group_sql,
                (composition_id, group_id_0))

            bucket = []
            for row in res:
                bucket.append((row[0],row[1]))

            iterator = permutations(bucket, r=2)

        else:

            res_0 = cur.execute(
                get_complex_group_sql,
                (composition_id, group_id_0))

            bucket_0 = []
            for row in res_0:
                bucket_0.append((row[0],row[1]))

            res_1 = cur.execute(
                get_complex_group_sql,
                (composition_id, group_id_1))

            bucket_1 = []
            for row in res_1:
                bucket_1.append((row[0],row[1]))

            iterator = product(bucket_0, bucket_1)



        for (reactants, products) in iterator:
            reaction = {
                'reactants' : reactants,
                'products' : products,
                'number_of_reactants' : len([i for i in reactants if i != -1]),
                'number_of_products' : len([i for i in products if i != -1])}


            decision_pathway = []
            if run_decision_tree(reaction,
                                 mol_entries,
                                 worker_payload.params,
                                 worker_payload.reaction_decision_tree,
                                 decision_pathway
                                 ):

                comm.send(
                    reaction,
                    dest=DISPATCHER_RANK,
                    tag=NEW_REACTION_DB)


            if run_decision_tree(reaction,
                                 mol_entries,
                                 worker_payload.params,
                                 worker_payload.logging_decision_tree):

                comm.send(
                    (reaction,
                     '\n'.join([str(f) for f in decision_pathway])
                     ),

                    dest=DISPATCHER_RANK,
                    tag=NEW_REACTION_LOGGING)
