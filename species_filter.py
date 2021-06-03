from mol_entry import MoleculeEntry


"""
Phase 1: species filtering
input: a list of dataset entries
output: a filtered list of mol_entries with fixed indices
description: this is where we remove isomorphic species, and do other forms of filtering. species_filter_decision_tree contains the filters
"""

def species_filter(dataset_entries, species_filter_decision_tree):
    mol_entries = [ MoleculeEntry.from_dataset_entry(e) for e in dataset_entries ]

    # TODO implement filtering

    for i, e in enumerate(mol_entries):
        e.parameters["ind"] = i

    return mol_entries
