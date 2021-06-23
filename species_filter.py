from mol_entry import MoleculeEntry
from functools import partial
from itertools import chain
from monty.serialization import dumpfn
import pickle
from species_questions import standard_mol_decision_tree, Terminal, run_decision_tree

"""
Phase 1: species filtering
input: a list of dataset entries
output: a filtered list of mol_entries with fixed indices
description: this is where we remove isomorphic species, and do other forms of filtering. Species decision tree is what we use for filtering.

species isomorphism filtering:

The input dataset entries will often contain isomorphic molecules. Identifying such isomorphisms doesn't fit into the species decision tree, so we have it as a preprocessing phase.
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
    l = map(partial(groupby, isomorphic_mols), chunks)

    return chain.from_iterable(l)


def species_filter(dataset_entries,
                   mol_entries_pickle_location,
                   species_decision_tree=standard_mol_decision_tree,
                   number_of_processes=8,
                   verbose=True):

    if verbose:
        print("starting species filter")
        print("loading molecule entries from json")

    mol_entries_unfiltered = [
        MoleculeEntry.from_dataset_entry(e) for e in dataset_entries ]


    # currently, take lowest energy mol in each iso class
    # if we want to add more non local species filtering it would go here

    if verbose:
        print("applying non local filters")

    def collapse_isomorphism_class(g):
        return min(g,key=lambda x: x.get_free_energy())


    mol_entries_no_iso = [
        collapse_isomorphism_class(g)
        for g in groupby_isomorphism(mol_entries_unfiltered, number_of_processes)]

    if verbose:
        print("applying local filters")

    # TODO: parallelize this
    mol_entries = [
        m for m in mol_entries_no_iso
        if run_decision_tree(m, species_decision_tree)]


    if verbose:
        print("assigning indices")

    for i, e in enumerate(mol_entries):
        e.parameters["ind"] = i


    if verbose:
        print("creating molecule entry pickle")
    # ideally we would serialize mol_entries to a json
    # some of the auxilary_data we compute
    # has frozen set keys, so doesn't seralize well into json format.
    # pickles work better in this setting
    with open(mol_entries_pickle_location, 'wb') as f:
        pickle.dump(mol_entries, f)


    return mol_entries
