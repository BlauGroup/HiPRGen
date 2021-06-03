from mol_entry import MoleculeEntry


"""
Phase 1: species filtering
input: a list of species
output: a filtered list of species with fixed indices
description: this is where we remove isomorphic species, and do other forms of filtering that we haven't thought of yet. We can organize things in such a way that adding more filters doesn't require someone to modify the inner loop
"""

def species_filter(dataset_entries):
    mol_entries = [ MoleculeEntry.from_dataset_entry(e) for e in dataset_entries ]

    # filters go here

    for i, e in enumerate(mol_entries):
        e.parameters["ind"] = i

    return mol_entries
