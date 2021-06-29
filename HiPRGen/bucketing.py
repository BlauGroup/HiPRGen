from mol_entry import MoleculeEntry
from itertools import combinations_with_replacement
import sqlite3

"""
Phase 2: bucketing pairs of species
input: filtered list of species with fixed indices
output: buckets labeled by atom count containing individual species and pairs of species
description: since each reaction conserves atom numbers, a concerted reaction only occurs between elements in a single bucket. There are tricks to reduce the number of pairs (like don't include (A,B) and (B,A)). If the number of species is 10,000, there are only 100 million such pairs which is within reach
"""

def create_table_sql(tag):
    return "CREATE TABLE " + tag + " (species_1, species_2)"

def insert_complex_sql(tag):
    return "INSERT INTO " + tag + " VALUES (?, ?)"

class Bucket:

    def insert_complex(self, tag, inds):
        cur = self.con.cursor()


        if tag not in self.insert_statements:
            cur.execute(create_table_sql(tag))
            self.insert_statements[tag] = insert_complex_sql(tag)

        cur.execute(self.insert_statements[tag], inds)

    def __init__(self, mol_entries, bucket_db, commit_freq=1000):
        self.con = sqlite3.connect(bucket_db)
        self.insert_statements = {}


        count = 0
        for m in mol_entries:
            tag = '_'.join(sorted(m.species))
            inds = (m.parameters["ind"], -1)
            self.insert_complex(tag, inds)

            count += 1
            if count % commit_freq == 0:
                self.con.commit()

        for (m1, m2) in combinations_with_replacement(mol_entries, 2):
            tag = '_'.join(sorted(m1.species + m2.species))
            inds = (m1.parameters["ind"], m2.parameters["ind"])
            self.insert_complex(tag, inds)

            count += 1
            if count % commit_freq == 0:
                self.con.commit()

        self.con.commit()
        self.con.close()

