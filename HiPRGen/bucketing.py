from HiPRGen.mol_entry import MoleculeEntry
from itertools import combinations_with_replacement
import sqlite3

"""
Phase 2: bucketing pairs of species input: filtered list of species
with fixed indices output: buckets labeled by atom count containing
individual species and pairs of species description: since each
reaction conserves atom numbers, a concerted reaction only occurs
between elements in a single bucket. There are tricks to reduce the
number of pairs (like don't include (A,B) and (B,A)). If the number of
species is 10,000, there are only 100 million such pairs which is
within reach
"""



def bucket(
        mol_entries,
        bucket_db,
        commit_freq=2000,
        group_size=1000):

    con = sqlite3.connect(bucket_db)
    cur = con.cursor()
    cur.execute(
        "CREATE TABLE complexes (species_1, species_2, composition_id, group_id)")

    # we create an index on (composition, group_id) so worker processes
    # during reaction filtering can read their work batch faster

    cur.execute(
        "CREATE INDEX composition_index ON complexes (composition_id, group_id)")

    group_counts = {}
    bucket_counts = {}
    composition_ids = {}
    commit_count = 0
    composition_count = 0

    for m in mol_entries:
        composition = '_'.join(sorted(m.species))

        if composition not in group_counts:
            group_counts[composition] = 0
            bucket_counts[composition] = 0
            composition_ids[composition] = composition_count
            composition_count += 1

        data = (m.ind, -1, composition_ids[composition], group_counts[composition])
        cur.execute("INSERT INTO complexes VALUES (?, ?, ?, ?)", data)

        commit_count += 1
        if commit_count % commit_freq == 0:
            con.commit()

        bucket_counts[composition] += 1
        if bucket_counts[composition] % group_size == 0:
            group_counts[composition] += 1


    for (m1, m2) in combinations_with_replacement(mol_entries, 2):
        composition = '_'.join(sorted(m1.species + m2.species))

        if composition not in group_counts:
            group_counts[composition] = 0
            bucket_counts[composition] = 0
            composition_ids[composition] = composition_count
            composition_count += 1


        data = (
            m1.ind,
            m2.ind,
            composition_ids[composition],
            group_counts[composition])

        cur.execute("INSERT INTO complexes VALUES (?, ?, ?, ?)", data)

        commit_count += 1
        if commit_count % commit_freq == 0:
            con.commit()

        bucket_counts[composition] += 1
        if bucket_counts[composition] % group_size == 0:
            group_counts[composition] += 1


    con.execute("CREATE TABLE group_counts (composition_id, count)")
    con.execute("CREATE TABLE compositions (composition_id, composition)")
    for composition in composition_ids:
        cur.execute(
            "INSERT INTO group_counts VALUES (?, ?)",
            (composition_ids[composition],
             group_counts[composition] + 1))

        cur.execute(
            "INSERT INTO compositions VALUES (?,?)",
            ((composition_ids[composition],
              composition)))



        commit_count += 1
        if commit_count % commit_freq == 0:
            con.commit()




    con.commit()
    con.close()


