from mol_entry import MoleculeEntry
from enum import Enum
from multiprocessing import Pool
from functools import partial
from itertools import chain


"""
Phase 1: species filtering
input: a list of dataset entries
output: a filtered list of mol_entries with fixed indices
description: this is where we remove isomorphic species, and do other forms of filtering. Species decision tree is what we use for filtering.

species isomorphism filtering:

The input dataset entries will often contain isomorphic molecules. Identifying such isomorphisms doesn't fit into the species decision tree, so we have it as a preprocessing phase.

species decision tree:

A question is a function q(mol_entry) -> Bool

Unlike for reaction filtering, these questions should not modify the mol_entry in any way.

A node is either a Terminal or a non empty list [(question, node)]

class Terminal(Enum):
    KEEP = 1
    DISCARD = -1

For the return value of a question, True means travel to this node and False means try next question in the list.

for non terminal nodes, it is an error if every question returns False. i.e getting stuck at a non terminal node is an error.

Once a Terminal node is reached, it tells us whether to keep or discard the species.

"""

def isomorphic_mols(m1, m2):
    return m1.mol_graph.isomorphic_to(m2.mol_graph)



def equal_tags(m1, m2):
    def tag(x):
        return ( x.formula, x.num_bonds, x.charge )

    return tag(m1) == tag(m2)


def groupby(equivalence_relation, xs):
    """
    warning: this has slightly different semantics than
    itertools groupby which depends on ordering.
    """
    groups = []

    for x in xs:
        group_found = False
        for group in groups:
            if equivalence_relation(x, group[0]):
                group.append(x)
                group_found = True
                break

        if not group_found:
            groups.append([x])

    return groups

def groupby_isomorphism(mol_entries, number_of_processes):
    """
    first group by tag to reduce the number of graph isos we need to compute
    """

    chunks = groupby(equal_tags, mol_entries)
    with Pool(number_of_processes) as p:
        l = p.map(partial(groupby, isomorphic_mols), chunks)

    return chain.from_iterable(l)



class Terminal(Enum):
    KEEP = 1
    DISCARD = -1

def run_decision_tree(mol_entry, decision_tree):
    node = decision_tree

    while type(node) == list:
        next_node = None
        for (question, new_node) in node:
            if question(mol_entry):
                next_node = new_node
                break

        node = next_node


    if type(node) == Terminal:
        if node == Terminal.KEEP:
            return True
        else:
            return False
    else:
        print(node)
        raise Exception("unexpected node type reached")


default_decision_tree = Terminal.KEEP

def species_filter(dataset_entries,
                   species_decision_tree,
                   number_of_processes=8):

    mol_entries_unfiltered = [
        MoleculeEntry.from_dataset_entry(e) for e in dataset_entries ]

    # currently, take lowest energy mol in each iso class
    def collapse_isomorphism_class(g):
        return min(g,key=lambda x: x.get_free_energy())

    mol_entries_no_iso = [
        collapse_isomorphism_class(g)
        for g in groupby_isomorphism(mol_entries_unfiltered, number_of_processes)]

    mol_entries = [
        m for m in mol_entries_no_iso
        if run_decision_tree(m, default_decision_tree)]


    for i, e in enumerate(mol_entries):
        e.parameters["ind"] = i

    return mol_entries
