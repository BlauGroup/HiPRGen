from HiPRGen.mol_entry import *
import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
import copy
from functools import partial
from HiPRGen.constants import *


"""
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

def run_decision_tree(mol_entry,
                      decision_tree,
                      decision_pathway=None):

    node = decision_tree

    while type(node) == list:
        next_node = None
        for (question, new_node) in node:
            if question(mol_entry):

                # if decision_pathway is a list,
                # append the question which
                # answered true i.e the edge we follow
                if decision_pathway is not None:
                    decision_pathway.append(question)


                next_node = new_node
                break

        node = next_node


    if type(node) == Terminal:
        if decision_pathway is not None:
            decision_pathway.append(node)


        if node == Terminal.KEEP:
            return True
        else:
            return False
    else:
        print(node)
        raise Exception("unexpected node type reached")



def metal_ion_filter(mol_entry):
    "only allow positively charged metal ions"
    if mol_entry.formula in m_formulas and mol_entry.charge <= 0:
        return True
    else:
        return False


def mol_not_connected(mol):
    return not nx.is_connected(mol.graph)


def add_star_hashes(mol):

    for i in range(mol.num_atoms):
        if i not in mol.m_inds:
            neighborhood = nx.generators.ego.ego_graph(
                mol.covalent_graph,
                i,
                1,
                undirected=True)

            mol.star_hashes[i] = weisfeiler_lehman_graph_hash(
                neighborhood,
                node_attr='specie')

    return False

def add_unbroken_fragment(width, mol):
    if mol.formula in m_formulas:
        return False

    neighborhood_hashes = {}

    for i in mol.non_metal_atoms:
        hash_list = []
        for d in range(width):
            neighborhood = nx.generators.ego.ego_graph(
                mol.covalent_graph,
                i,
                d,
                undirected=True)

            neighborhood_hash = weisfeiler_lehman_graph_hash(
                neighborhood,
                node_attr='specie')

            hash_list.append(neighborhood_hash)

        neighborhood_hashes[i] = hash(tuple(hash_list))

    fragment = Fragment(
        mol.covalent_hash,
        mol.non_metal_atoms,
        neighborhood_hashes,
        mol.covalent_graph,
        []
    )

    fragment_complex = FragmentComplex(
        1,
        0,
        [],
        [fragment])

    mol.fragment_data.append(fragment_complex)

    return False


def add_single_bond_fragments(width, mol):

    if mol.formula in m_formulas:
        return False



    for edge in mol.covalent_graph.edges:
        fragments = []
        h = copy.deepcopy(mol.covalent_graph)
        h.remove_edge(*edge)
        connected_components = nx.algorithms.components.connected_components(h)
        for c in connected_components:

            subgraph = h.subgraph(c)

            fragment_hash = weisfeiler_lehman_graph_hash(
                subgraph,
                node_attr='specie')


            neighborhood_hashes = {}
            for i in c:
                hash_list = []
                for d in range(width):
                    neighborhood = nx.generators.ego.ego_graph(
                        subgraph,
                        i,
                        d,
                        undirected=True)

                    neighborhood_hash = weisfeiler_lehman_graph_hash(
                        neighborhood,
                        node_attr='specie')

                    hash_list.append(neighborhood_hash)

                neighborhood_hashes[i] = hash(tuple(hash_list))

            fragment = Fragment(
                fragment_hash,
                c,
                neighborhood_hashes,
                subgraph,
                [i for i in c if i in edge[0:2]]
            )

            fragments.append(fragment)

        fragment_complex = FragmentComplex(
            len(fragments),
            1,
            [edge[0:2]],
            fragments)

        mol.fragment_data.append(fragment_complex)

    return False


def metal_complex(mol):
    # if mol is a metal, it isn't a metal complex
    if mol.formula in m_formulas:
        return False

    return not nx.is_connected(mol.covalent_graph)


def li_fix_hydrogen_bonding(mol):


    if mol.num_atoms > 1:
        for i in range(mol.num_atoms):
            if mol.species[i] == 'H':

                adjacent_atoms = []

                for bond in mol.graph.edges:
                    if i in bond[0:2]:

                        if i == bond[0]:
                            adjacent_atom = bond[1]
                        else:
                            adjacent_atom = bond[0]

                        displacement = (mol.atom_locations[adjacent_atom] -
                                        mol.atom_locations[i])

                        dist = np.inner(displacement, displacement)

                        adjacent_atoms.append((adjacent_atom, dist))


                closest_atom, _ = min(adjacent_atoms, key=lambda pair: pair[1])

                for adjacent_atom, _ in adjacent_atoms:
                    if adjacent_atom != closest_atom:
                        mol.graph.remove_edge(i, adjacent_atom)
                        mol.covalent_graph.remove_edge(i, adjacent_atom)

    return False


def mg_fix_hydrogen_bonding(mol):

    # TODO: Make this user-defined?
    max_dist = 1.5

    if mol.num_atoms > 1:
        for i in range(mol.num_atoms):
            if mol.species[i] == 'H':

                adjacent_atoms = []

                for bond in mol.graph.edges:
                    if i in bond[0:2]:

                        if i == bond[0]:
                            adjacent_atom = bond[1]
                        else:
                            adjacent_atom = bond[0]

                        displacement = (mol.atom_locations[adjacent_atom] -
                                        mol.atom_locations[i])

                        dist = np.inner(displacement, displacement)

                        adjacent_atoms.append((adjacent_atom, dist))

                for adjacent_atom, dist in adjacent_atoms:
                    if dist > max_dist ** 2:
                        mol.graph.remove_edge(i, adjacent_atom)
                        if adjacent_atom in mol.covalent_graph:
                            mol.covalent_graph.remove_edge(i, adjacent_atom)

    return False


def bad_metal_coordination(mol):

    if mol.formula not in m_formulas:

        if (len(metals.intersection(set(mol.species))) > 0 and
            mol.number_of_coordination_bonds == 0):

            return True

    return False

def li_set_solvation_free_energy(
        solvation_env,
        mol):
    """
    metal atoms coordinate with the surrounding solvent. We need to correct
    free energy to take this into account. The correction is
    solvation_correction * (
           max_coodination_bonds -
           number_of_coordination_bonds_in_mol).
    Since coordination bonding can't reliably be detected from the molecule
    graph, we search for all atoms within a radius of the metal atom and
    discard them if they are positively charged.
    """

    mol.number_of_coordination_bonds = 0

    correction = 0.0

    for i in mol.m_inds:

        species = mol.species[i]
        coordination_partners = []
        radius = solvation_env["coordination_radius"][species]

        for j in range(mol.num_atoms):
            if j != i:
                displacement_vector = (
                    mol.atom_locations[j] -
                    mol.atom_locations[i])

                if (np.inner(displacement_vector, displacement_vector)
                    < radius ** 2 and (
                        mol.partial_charges_resp[j] < 0 or
                        mol.partial_charges_mulliken[j] < 0)):
                    coordination_partners.append(j)


        number_of_coordination_bonds = len(coordination_partners)
        mol.number_of_coordination_bonds += number_of_coordination_bonds
        correction += solvation_env["solvation_correction"][species] * (
            solvation_env["max_number_of_coordination_bonds"][species] -
            number_of_coordination_bonds)

    mol.solvation_free_energy = correction + mol.free_energy
    return False


def mg_set_solvation_free_energy(
        solvation_env,
        mol):
    """
    metal atoms coordinate with the surrounding solvent. We need to correct
    free energy to take this into account. The correction is
    solvation_correction * (
           max_coodination_bonds -
           number_of_coordination_bonds_in_mol).
    Since coordination bonding can't reliably be detected from the molecule
    graph, we search for all atoms within a radius of the metal atom and
    discard them if they are positively charged.
    """

    correction = 0.0
    mol.number_of_coordination_bonds = 0

    for i in mol.m_inds:

        species = mol.species[i]
        partial_charge = mol.partial_charges_mulliken[i]

        if partial_charge < 1.2:
            effective_charge = "_1"
        elif partial_charge >= 1.2:
            effective_charge = "_2"


        coordination_partners = list()
        species_charge = species + effective_charge
        radius = solvation_env["coordination_radius"][species_charge]

        for j in range(mol.num_atoms):
            if j != i:
                displacement_vector = (
                    mol.atom_locations[j] -
                    mol.atom_locations[i])
                if (np.inner(displacement_vector, displacement_vector)
                    < radius ** 2 and (
                        mol.partial_charges_resp[j] < 0 or
                        mol.partial_charges_mulliken[j] < 0)):
                    coordination_partners.append(j)

        number_of_coordination_bonds = len(coordination_partners)
        mol.number_of_coordination_bonds += number_of_coordination_bonds
        correction += solvation_env["solvation_correction"][species_charge] * (
            solvation_env["max_number_of_coordination_bonds"][species_charge] -
            number_of_coordination_bonds)

    mol.solvation_free_energy =  correction + mol.free_energy
    return False



def default_true(mol):
    return True


fragment_neighborhood_width = 5

# any species filter which modifies bonding has to come before
# any filter checking for connectivity (which includes the metal-centric complex filter)

li_species_decision_tree = [
    (partial(li_set_solvation_free_energy,
             li_ec), Terminal.KEEP),

    (li_fix_hydrogen_bonding, Terminal.KEEP),
    (metal_ion_filter, Terminal.DISCARD),
    (bad_metal_coordination, Terminal.DISCARD),
    (mol_not_connected, Terminal.DISCARD),
    (metal_complex, Terminal.DISCARD),
    (add_star_hashes, Terminal.KEEP),

    (partial(add_unbroken_fragment,
             fragment_neighborhood_width), Terminal.KEEP),
    (partial(add_single_bond_fragments,
             fragment_neighborhood_width), Terminal.KEEP),

    (default_true, Terminal.KEEP)
    ]

mg_thf_species_decision_tree = [
    (partial(mg_set_solvation_free_energy,
             mg_thf), Terminal.KEEP),

    (mg_fix_hydrogen_bonding, Terminal.KEEP),
    (metal_ion_filter, Terminal.DISCARD),
    (bad_metal_coordination, Terminal.DISCARD),
    (mol_not_connected, Terminal.DISCARD),
    (metal_complex, Terminal.DISCARD),
    (add_star_hashes, Terminal.KEEP),

    (partial(add_unbroken_fragment,
             fragment_neighborhood_width), Terminal.KEEP),
    (partial(add_single_bond_fragments,
             fragment_neighborhood_width), Terminal.KEEP),

    (default_true, Terminal.KEEP)
    ]

standard_species_logging_decision_tree = Terminal.KEEP
