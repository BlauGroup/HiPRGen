from mol_entry import MoleculeEntry


"""
Phase 1: species filtering
input: a list of dataset entries
output: a filtered list of mol_entries with fixed indices
description: this is where we remove isomorphic species, and do other forms of filtering. Species decision tree is what we use for filtering.

species decision tree

A question is a function q(mol_entry) -> Bool

A node is either a Terminal or a non empty list [(question, node)]

class Terminal(Enum):
    KEEP = 1
    DISCARD = -1

For the return value of a question, True means travel to this node and False means try next question in the list.

for non terminal nodes, it is an error if every question returns False. i.e getting stuck at a non terminal node is an error.

Once a Terminal node is reached, it tells us whether to keep or discard the species.

"""

def species_filter(dataset_entries, species_decision_tree):
    mol_entries = [ MoleculeEntry.from_dataset_entry(e) for e in dataset_entries ]

    # TODO implement filtering

    for i, e in enumerate(mol_entries):
        e.parameters["ind"] = i

    return mol_entries
