from mol_entry import MoleculeEntry
from itertools import combinations_with_replacement

"""
Phase 2: bucketing pairs of species
input: filtered list of species with fixed indices
output: buckets labeled by atom count containing individual species and pairs of species
description: since each reaction conserves atom numbers, a concerted reaction only occurs between elements in a single bucket. There are tricks to reduce the number of pairs (like don't include (A,B) and (B,A)). If the number of species is 10,000, there are only 100 million such pairs which is within reach
"""

def bucket(mol_entries, bucket_folder):
    for (m1, m2) in combinations_with_replacement(mol_entries, 2):
        tag = '_'.join(sorted(m1.species + m2.species))
        print(tag)
