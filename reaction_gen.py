from mol_entry import MoleculeEntry
from itertools import combinations
from multiprocessing import Process, Queue
import sqlite3

"""
Phase 3: reaction gen and filtering
input: a bucket labeled by atom count
output: a list of reactions from that bucket
description: Loop through all possible reactions in the bucket and apply a list of filters. This will run in parallel over each bucket. The filters are organized into a decision tree so we don't need to run expensive filters for every single reaction (the expensive ones occur deep in the tree). Filters will be writable by hand, or we could have machine learning filters. Also, the reaction rate will be set during this phase so need to percolate through the decision tree

Phase 4: collating and indexing
input: all the outputs of phase 3 as they are generated (so phase 3 and 4 will be running at the same time)
output: the final list of reactions (or db assuming it is to big to hold in memory)
description: the worker processes from phase 3 are sending their reactions to this phase and it is writing them to DB as it gets them. We can ensure that duplicates don't get generated in phase 3, so we will actually get a massive performance boost for writing to the DB also
"""

def list_or(a_list):
    return True in a_list

create_metadata_table = """
    CREATE TABLE metadata (
            number_of_species   INTEGER NOT NULL,
            number_of_reactions INTEGER NOT NULL
    );
"""

insert_metadata = """
    INSERT INTO metadata VALUES (?, ?)
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

def dispatcher(mol_entries, bucket_db, rn_db, commit_freq=1000):
    reaction_queue = Queue()

    bucket_con = sqlite3.connect(bucket_db)
    bucket_cur = bucket_con.cursor()

    processes = {}
    res = bucket_cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    for name in res:
        table = name[0]
        p = Process(
            target=reaction_filter,
            args=(
                mol_entries,
                bucket_db,
                table,
                reaction_queue))

        processes[table] = p

    rn_con = sqlite3.connect(rn_db)
    rn_cur = rn_con.cursor()
    rn_cur.execute(create_metadata_table)
    rn_cur.execute(create_reactions_table)
    rn_con.commit()

    for table in processes:
        processes[table].start()

    living_children = True
    reaction_index = 0

    while living_children:
        if reaction_queue.empty():
            living_bools = [processes[table].is_alive() for table in processes]
            living_children = list_or(living_bools)

        else:
            reaction = reaction_queue.get()
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

    rn_cur.execute(
        insert_metadata,
        (len(mol_entries) + 1,
         reaction_index + 1))

    rn_con.commit()
    bucket_con.close()
    rn_con.close()




def reaction_filter(mol_entries, bucket_db, table, reaction_queue):
    con = sqlite3.connect(bucket_db)
    cur = con.cursor()
    bucket = []


    res = cur.execute("SELECT * FROM " + table)
    for pair in res:
        bucket.append(pair)

    for (reactants, products) in combinations(bucket, 2):
        reaction = {
            'reactants' : reactants,
            'products' : products,
            'number_of_reactants' : len([i for i in reactants if i != -1]),
            'number_of_products' : len([i for i in products if i != -1]),
            'rate' : 0.0,
            'dG' : 0.0 }



        reaction_queue.put(reaction)

