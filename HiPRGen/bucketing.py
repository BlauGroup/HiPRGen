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
    composition_count = 0

    # Rows are buffered and inserted with executemany rather than one
    # execute() per row, since for N species the two loops below insert
    # O(N^2) rows and per-row execute() calls make that count of
    # Python<->sqlite round trips the dominant cost of this phase.
    # Insertion order doesn't affect the table's contents, so batching
    # doesn't change what ends up in the database, only how it gets there.
    insert_buffer = []

    def flush():
        if insert_buffer:
            cur.executemany(
                "INSERT INTO complexes VALUES (?, ?, ?, ?)", insert_buffer)
            con.commit()
            insert_buffer.clear()

    # Each molecule's species list is sorted exactly once up front, rather
    # than resorting the same molecule's species list from scratch on
    # every one of the ~N pairings it appears in below. The composition
    # key is still built as a '_'.join()-ed string, not a tuple: a
    # composition's dict entries get touched 2-3 times per pair below
    # (membership check, bucket_counts increment, sometimes group_counts
    # increment), and CPython caches a string's hash on the object after
    # the first of those touches, making the rest free - a tuple has no
    # such cache and recomputes its combined hash on every touch, which
    # measured out to a net wash despite being cheaper to construct.
    sorted_species = [sorted(m.species) for m in mol_entries]
    paired = list(zip(mol_entries, sorted_species))

    for m, ss in paired:
        composition = '_'.join(ss)

        if composition not in group_counts:
            group_counts[composition] = 0
            bucket_counts[composition] = 0
            composition_ids[composition] = composition_count
            composition_count += 1

        insert_buffer.append(
            (m.ind, -1, composition_ids[composition], group_counts[composition]))
        if len(insert_buffer) >= commit_freq:
            flush()

        bucket_counts[composition] += 1
        if bucket_counts[composition] % group_size == 0:
            group_counts[composition] += 1


    for ((m1, ss1), (m2, ss2)) in combinations_with_replacement(paired, 2):
        composition = '_'.join(sorted(ss1 + ss2))

        if composition not in group_counts:
            group_counts[composition] = 0
            bucket_counts[composition] = 0
            composition_ids[composition] = composition_count
            composition_count += 1

        insert_buffer.append((
            m1.ind,
            m2.ind,
            composition_ids[composition],
            group_counts[composition]))
        if len(insert_buffer) >= commit_freq:
            flush()

        bucket_counts[composition] += 1
        if bucket_counts[composition] % group_size == 0:
            group_counts[composition] += 1

    flush()

    con.execute("CREATE TABLE group_counts (composition_id, count)")
    con.execute("CREATE TABLE compositions (composition_id, composition)")

    cur.executemany(
        "INSERT INTO group_counts VALUES (?, ?)",
        [
            (composition_ids[composition], group_counts[composition] + 1)
            for composition in composition_ids
        ])

    cur.executemany(
        "INSERT INTO compositions VALUES (?, ?)",
        [
            (composition_ids[composition], composition)
            for composition in composition_ids
        ])

    con.commit()
    con.close()


