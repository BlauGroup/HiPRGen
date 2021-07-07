from HiPRGen.mol_entry import MoleculeEntry
from functools import partial
from itertools import chain
from monty.serialization import dumpfn
import pickle
from HiPRGen.species_questions import standard_mol_decision_tree, Terminal, run_decision_tree
from time import localtime, strftime
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash

"""
Phase 1: species filtering
input: a list of dataset entries
output: a filtered list of mol_entries with fixed indices
description: this is where we remove isomorphic species, and do other forms of filtering. Species decision tree is what we use for filtering.

species isomorphism filtering:

The input dataset entries will often contain isomorphic molecules. Identifying such isomorphisms doesn't fit into the species decision tree, so we have it as a preprocessing phase.
"""

def groupby_isomorphism(mols):
    isomorphism_buckets = {}
    for mol in mols:
        mol_hash = weisfeiler_lehman_graph_hash(
            mol.graph.to_undirected(),
            node_attr='specie'
        )

        tag = (mol.charge, mol.formula, mol.num_bonds, mol_hash)

        if tag in isomorphism_buckets:
            isomorphism_buckets[tag].append(mol)
        else:
            isomorphism_buckets[tag] = [mol]

    return isomorphism_buckets

def log_message(string):
    print(
        '[' + strftime('%H:%M:%S', localtime()) + ']',
        string)

def species_filter(dataset_entries,
                   mol_entries_pickle_location,
                   species_decision_tree=standard_mol_decision_tree
                   ):

    log_message("starting species filter")
    log_message("loading molecule entries from json")

    mol_entries_unfiltered = [
        MoleculeEntry.from_dataset_entry(e) for e in dataset_entries ]


    # currently, take lowest energy mol in each iso class
    # if we want to add more non local species filtering it would go here

    log_message("applying non local filters")

    def collapse_isomorphism_class(g):
        return min(g,key=lambda x: x.get_free_energy())


    mol_entries_no_iso = [
        collapse_isomorphism_class(g)
        for g in groupby_isomorphism(mol_entries_unfiltered).values()]

    log_message("applying local filters")

    mol_entries = [
        m for m in mol_entries_no_iso
        if run_decision_tree(m, species_decision_tree)]


    log_message("assigning indices")

    for i, e in enumerate(mol_entries):
        e.parameters["ind"] = i


    log_message("creating molecule entry pickle")
    # ideally we would serialize mol_entries to a json
    # some of the auxilary_data we compute
    # has frozen set keys, so doesn't seralize well into json format.
    # pickles work better in this setting
    with open(mol_entries_pickle_location, 'wb') as f:
        pickle.dump(mol_entries, f)

    return mol_entries
