from HiPRGen.mol_entry import MoleculeEntry
from functools import partial
from itertools import chain
from monty.serialization import dumpfn
import pickle
from HiPRGen.species_questions import (
    run_decision_tree,
    standard_species_logging_decision_tree
)
import networkx as nx
from time import localtime, strftime
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
import networkx.algorithms.isomorphism as iso
from HiPRGen.report_generator import ReportGenerator
"""
Phase 1: species filtering
input: a list of dataset entries
output: a filtered list of mol_entries with fixed indices
description: this is where we remove isomorphic species, and do other forms of filtering. Species decision tree is what we use for filtering.

species isomorphism filtering:

The input dataset entries will often contain isomorphic molecules. Identifying such isomorphisms doesn't fit into the species decision tree, so we have it as a preprocessing phase.
"""

def sort_into_tags(mols):
    isomorphism_buckets = {}
    for mol in mols:

        tag = (mol.charge, mol.formula, mol.covalent_hash)

        if tag in isomorphism_buckets:
            isomorphism_buckets[tag].append(mol)
        else:
            isomorphism_buckets[tag] = [mol]

    return isomorphism_buckets


def really_covalent_isomorphic(mol1, mol2):
    """
    check for isomorphism directly instead of using hash.
    warning: this is really slow. It is used in species filtering
    to avoid hash collisions. Do not use it anywhere else.
    """
    return nx.is_isomorphic(
        mol1.covalent_graph,
        mol2.covalent_graph,
        node_match = iso.categorical_node_match('specie', None)
    )



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


def log_message(string):
    print(
        '[' + strftime('%H:%M:%S', localtime()) + ']',
        string)

def species_filter(
        dataset_entries,
        mol_entries_pickle_location,
        species_report,
        species_decision_tree,
        coordimer_weight,
        species_logging_decision_tree=standard_species_logging_decision_tree,
        generate_unfiltered_mol_pictures=False):

    """
    run each molecule through the species decision tree and then choose the lowest weight
    coordimer based on the coordimer_weight function.
    """

    log_message("starting species filter")
    log_message("loading molecule entries from json")

    mol_entries_unfiltered = [
        MoleculeEntry.from_dataset_entry(e) for e in dataset_entries ]


    log_message("generating unfiltered mol pictures")

    report_generator = ReportGenerator(
        mol_entries_unfiltered,
        species_report,
        mol_pictures_folder_name='mol_pictures_unfiltered',
        rebuild_mol_pictures=generate_unfiltered_mol_pictures
    )

    report_generator.emit_text("species report")

    log_message("applying local filters")
    mol_entries_filtered = []

    # note: it is important here that we are applying the local filters before
    # the non local ones. We remove some molecules which are lower energy
    # than other more realistic lithomers.

    for i, mol in enumerate(mol_entries_unfiltered):
        log_message("filtering " + mol.entry_id)
        decision_pathway = []
        if run_decision_tree(mol, species_decision_tree, decision_pathway):
            mol_entries_filtered.append(mol)

        if run_decision_tree(mol, species_logging_decision_tree):

            report_generator.emit_verbatim(
                '\n'.join([str(f) for f in decision_pathway]))


            report_generator.emit_text("entry id: " + mol.entry_id)
            report_generator.emit_text("uncorrected free energy: " +
                                       str(mol.free_energy))

            report_generator.emit_text(
                "number of coordination bonds: " +
                str(mol.number_of_coordination_bonds))

            report_generator.emit_text(
                "corrected free energy: " +
                str(mol.solvation_free_energy))

            report_generator.emit_text(
                "formula: " + mol.formula)

            report_generator.emit_molecule(i, include_index=False)
            report_generator.emit_newline()


    report_generator.finished()


    # python doesn't have shared memory. That means that every worker during
    # reaction filtering must maintain its own copy of the molecules.
    # for this reason, it is good to remove attributes that are only used
    # during species filtering.
    log_message("clearing unneeded attributes")
    for m in mol_entries_filtered:
        m.partial_charges_resp = None
        m.partial_charges_mulliken = None
        m.partial_charges_nbo = None
        m.atom_locations = None

    # currently, take lowest energy mol in each iso class
    log_message("applying non local filters")

    # when we choose a coordimer, we also keep track of the others
    # so we can use them to compute redox rates
    def collapse_isomorphism_group(g):
        lowest_energy_coordimer = min(g,key=coordimer_weight)


        if len(lowest_energy_coordimer.m_inds) > 0:
            coordimers = {}

            for m in g:
                coordimers[m.total_hash] = m

            lowest_energy_coordimer.coordimers = coordimers

        return lowest_energy_coordimer


    mol_entries = []

    for tag_group in sort_into_tags(mol_entries_filtered).values():
        for iso_group in groupby(really_covalent_isomorphic, tag_group):
            mol_entries.append(
                collapse_isomorphism_group(iso_group))


    log_message("assigning indices")

    for i, e in enumerate(mol_entries):
        e.ind = i


    log_message("creating molecule entry pickle")
    # ideally we would serialize mol_entries to a json
    # some of the auxilary_data we compute
    # has frozen set keys, so doesn't seralize well into json format.
    # pickles work better in this setting
    with open(mol_entries_pickle_location, 'wb') as f:
        pickle.dump(mol_entries, f)

    log_message("species filtering finished. " +
                str(len(mol_entries)) +
                " species")

    return mol_entries
