from HiPRGen.mol_entry import MoleculeEntry, find_fragment_atom_mappings, build_compressed_graph
from functools import partial
from itertools import chain
from monty.serialization import dumpfn
import pickle
import copy
from HiPRGen.species_questions import run_decision_tree
from HiPRGen.constants import Terminal
import networkx as nx
from time import localtime, strftime
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
import networkx.algorithms.isomorphism as iso
from HiPRGen.report_generator import ReportGenerator
from pymatgen.core.periodic_table import DummySpecies
from pymatgen.core.sites import Site
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from bondnet.model.training_utils import get_grapher
from bondnet.core.molwrapper import MoleculeWrapper
from bondnet.data.transformers import HeteroGraphFeatureStandardScaler


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
        node_match=iso.categorical_node_match("specie", None),
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
    print("[" + strftime("%H:%M:%S", localtime()) + "]", string)


def species_filter(
    dataset_entries,
    mol_entries_pickle_location,
    dgl_mol_grphs_pickle_location,
    grapher_features_pickle_location,
    species_report,
    species_decision_tree,
    coordimer_weight,
    species_logging_decision_tree=Terminal.DISCARD,
    generate_unfiltered_mol_pictures=False,
):

    """
    run each molecule through the species decision tree and then choose the lowest weight
    coordimer based on the coordimer_weight function.
    """

    log_message("starting species filter")
    log_message("loading molecule entries from json")

    if "has_props" in dataset_entries[0].keys():
        mol_entries_unfiltered = [MoleculeEntry.from_mp_doc(e) for e in dataset_entries]
        log_message("MP doc entries passed")
    else:
        log_message("dataset entries passed")
        mol_entries_unfiltered = [
            MoleculeEntry.from_dataset_entry(e) for e in dataset_entries
        ]

    log_message("found " + str(len(mol_entries_unfiltered)) + " molecule entries")
    log_message("generating unfiltered mol pictures")

    report_generator = ReportGenerator(
        mol_entries_unfiltered,
        species_report,
        mol_pictures_folder_name="mol_pictures_unfiltered",
        rebuild_mol_pictures=generate_unfiltered_mol_pictures,
    )

    report_generator.emit_text("species report")

    log_message("applying local filters")
    mol_entries_filtered = []
    elements = set()
    # note: it is important here that we are applying the local filters before
    # the non local ones. We remove some molecules which are lower energy
    # than other more realistic lithomers.

    for i, mol in enumerate(mol_entries_unfiltered):
        log_message("filtering " + mol.entry_id)
        decision_pathway = []
        if run_decision_tree(mol, species_decision_tree, decision_pathway):
            mol_entries_filtered.append(mol)

        if run_decision_tree(mol, species_logging_decision_tree):
            # create a set of elements
            # log_message("create a set of elements")
            for element in mol.species:
                if element not in elements:
                    elements.add(element)

            report_generator.emit_verbatim(
                "\n".join([str(f) for f in decision_pathway])
            )

            report_generator.emit_text("number: " + str(i))
            report_generator.emit_text("entry id: " + mol.entry_id)
            report_generator.emit_text(
                "uncorrected free energy: " + str(mol.free_energy)
            )

            # report_generator.emit_text(
            #     "number of coordination bonds: " + str(mol.number_of_coordination_bonds)
            # )

            # report_generator.emit_text(
            #     "corrected free energy: " + str(mol.solvation_free_energy)
            # )

            report_generator.emit_text("formula: " + mol.formula)

            report_generator.emit_molecule(i, include_index=False)
            report_generator.emit_newline()

    report_generator.finished()

    # python doesn't have shared memory. That means that every worker during
    # reaction filtering must maintain its own copy of the molecules.
    # for this reason, it is good to remove attributes that are only used
    # during species filtering.
    log_message("clearing unneeded attributes")
    for m in mol_entries_filtered:
        del m.partial_charges_resp
        del m.partial_charges_mulliken
        del m.partial_charges_nbo
        del m.partial_spins_nbo
        del m.atom_locations

    # currently, take lowest energy mol in each iso class
    log_message("applying non local filters")

    def collapse_isomorphism_group(g):
        lowest_energy_coordimer = min(g, key=coordimer_weight)
        return lowest_energy_coordimer

    mol_entries = []

    for tag_group in sort_into_tags(mol_entries_filtered).values():
        for iso_group in groupby(really_covalent_isomorphic, tag_group):
            mol_entries.append(collapse_isomorphism_group(iso_group))

    log_message("assigning indices")

    for i, e in enumerate(mol_entries):
        e.ind = i


    log_message("mapping fragments")
    fragment_dict = {}
    for mol in mol_entries:
        # print(mol.entry_id)
        for fragment_complex in mol.fragment_data:
            for ii, fragment in enumerate(fragment_complex.fragment_objects):
                hot_nbh_hashes = list(fragment.hot_atoms.keys())
                assert len(hot_nbh_hashes) == 0 or len(hot_nbh_hashes) == 1 or len(hot_nbh_hashes) == 2
                assert fragment.fragment_hash == fragment_complex.fragment_hashes[ii]
                if fragment.fragment_hash not in fragment_dict:
                    fragment_dict[fragment.fragment_hash] = copy.deepcopy(fragment)
                all_mappings = find_fragment_atom_mappings(
                    fragment,
                    fragment_dict[fragment.fragment_hash])
                fragment_complex.fragment_mappings.append(all_mappings)
        assert len(fragment_complex.fragment_objects) == len(fragment_complex.fragment_mappings)


    log_message(str(len(fragment_dict.keys())) + " unique fragments found")

    # Make DGL Molecule graphs via BonDNet functions
    log_message("creating dgl molecule graphs")
    dgl_molecules_dict = {}
    dgl_molecules = []
    extra_keys = []
    
    for mol in mol_entries:
        # print(f"mol: {mol.mol_graph}")
        molecule_grapher = get_grapher(extra_keys)

        non_metal_bonds = [ tuple(sorted([i, j])) for i, j, _ in mol.covalent_graph.edges.data()]
        # print(f"non metal bonds: {non_metal_bonds}")
        mol_wrapper = MoleculeWrapper(mol_graph = mol.mol_graph, free_energy = mol.energy, id = mol.entry_id, non_metal_bonds = non_metal_bonds)
        feature = {'charge': mol.charge}
        dgl_molecule_graph = molecule_grapher.build_graph_and_featurize(mol_wrapper, extra_feats_info = feature, dataset_species = elements)
        dgl_molecules.append(dgl_molecule_graph)
        dgl_molecules_dict[mol.entry_id] = mol.ind
    grapher_features= {'feature_size':molecule_grapher.feature_size, 'feature_name': molecule_grapher.feature_name}
    #mol_wrapper_dict[mol.entry_id] = mol_wrapper

    # Normalize DGL molecule graphs
    scaler = HeteroGraphFeatureStandardScaler(mean = None, std = None)
    normalized_graphs = scaler(dgl_molecules)

    # print(f"mean: {scaler._mean}")
    # print(f"std: {scaler._std}")
    

    # Create a dictionary where key is mol.entry_id and value is a normalized dgl molecule graph
    for key in dgl_molecules_dict.keys():
        temp_index = dgl_molecules_dict[key]
        dgl_molecules_dict[key] = normalized_graphs[temp_index]
    #print(dgl_molecules_dict)


    log_message("creating molecule entry pickle")
    # ideally we would serialize mol_entries to a json
    # some of the auxilary_data we compute
    # has frozen set keys, so doesn't seralize well into json format.
    # pickles work better in this setting
    with open(mol_entries_pickle_location, "wb") as f:
        pickle.dump(mol_entries, f)

    with open(dgl_mol_grphs_pickle_location, "wb") as f:
        pickle.dump(dgl_molecules_dict, f)
    
    with open(grapher_features_pickle_location, "wb") as f:
        pickle.dump(grapher_features, f)

    log_message("species filtering finished. " + str(len(mol_entries)) + " species")

    return mol_entries, dgl_molecules_dict


def add_electron_species(
    mol_entries, mol_entries_pickle_location, electron_free_energy
):
    e_site = Site(
        DummySpecies("E", oxidation_state=None, properties=None), [0.0, 0.0, 0.0]
    )
    e_mol = Molecule.from_sites([e_site])
    e_mol.set_charge_and_spin(-1, 2)
    e_graph = MoleculeGraph.with_empty_graph(molecule=e_mol)
    electron_entry = MoleculeEntry(
        molecule=e_mol,
        energy=electron_free_energy,
        enthalpy=0,
        entropy=0,
        entry_id=None,
        mol_graph=e_graph,
        partial_charges_resp=None,
        partial_charges_mulliken=None,
        partial_charges_nbo=None,
        electron_affinity=None,
        ionization_energy=None,
        spin_multiplicity=None,
        partial_spins_nbo=None,
    )
    electron_entry.ind = len(mol_entries)
    mol_entries.append(electron_entry)
    with open(mol_entries_pickle_location, "wb") as f:
        pickle.dump(mol_entries, f)
    return mol_entries
