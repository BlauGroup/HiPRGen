from HiPRGen.mol_entry import MoleculeEntry, FragmentComplex, FragmentObject
import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
import copy
from functools import partial
from HiPRGen.constants import li_ec, Terminal, mg_g2, mg_thf, m_formulas, metals
import numpy as np
from monty.json import MSONable
from itertools import combinations

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


def run_decision_tree(mol_entry, decision_tree, decision_pathway=None):

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


class metal_ion_filter(MSONable):
    "only allow positively charged metal ions"

    def __init__(self):
        pass

    def __call__(self, mol_entry):
        if mol_entry.formula in m_formulas and mol_entry.charge <= 0:
            return True
        else:
            return False


class mol_not_connected(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        return not nx.is_connected(mol.graph)


class spin_multiplicity_filter(MSONable):
    def __init__(self, threshold):
        self.threshold = threshold

    def __call__(self, mol):
        if mol.spin_multiplicity == 2:
            num_partial_spins_above_threshold = 0
            for i in range(mol.num_atoms):
                if mol.partial_spins_nbo[i] > self.threshold:
                    num_partial_spins_above_threshold += 1

            if num_partial_spins_above_threshold >= 2:
                mol.penalty += 1

        return False


class positive_penalty(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        if mol.penalty > 0:
            return True
        else:
            return False


class add_star_hashes(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        for i in range(mol.num_atoms): #iterates over all atoms in a molecule
            if i not in mol.m_inds: #ignoring metal atoms
                neighborhood = nx.generators.ego.ego_graph(   #generates an ego graph named neighborhood, with atom i at the
                    mol.covalent_graph, i, 1, undirected=True #center, with the nodes of the graph being the atoms i is
                )                                             #covalently bonded to

                mol.star_hashes[i] = weisfeiler_lehman_graph_hash( #star_hashes is a dictionary, and this adds an entry to it
                    neighborhood, node_attr="specie"               #with the atom index, i, as the key and a graph_hash (string)
                )                                                  #as the value

        return False


class add_unbroken_fragment(MSONable):  #aka adds unfragmented molecule as a "fragment complex"
    def __init__(self, neighborhood_width=None):
        self.neighborhood_width = neighborhood_width

    def __call__(self, mol):
        if mol.formula in m_formulas:
            return False

        fragment_objects = []
        if self.neighborhood_width is not None:
            neighborhood_hashes = {}
            for i in mol.non_metal_atoms:
                hash_list = []
                for d in range(self.neighborhood_width):
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
            fragment_object = FragmentObject(
                mol.covalent_hash,
                mol.non_metal_atoms,
                neighborhood_hashes,
                mol.covalent_graph,
                []
            )
            fragment_objects.append(fragment_object)

        fragment_complex = FragmentComplex(1, 0, [], [mol.covalent_hash], fragment_objects)

        mol.fragment_data.append(fragment_complex)

        return False


class add_single_bond_fragments(MSONable): #called for all species that have passed through filtration
    def __init__(self, allow_ring_opening=True, neighborhood_width=None):
        self.allow_ring_opening = allow_ring_opening
        self.neighborhood_width = neighborhood_width

    def __call__(self, mol):

        if mol.formula in m_formulas:
            return False

        for edge in mol.covalent_graph.edges: #iterates through each bond in a molecule graph by iterating through a list of tuples
            fragment_hashes = []
            fragment_objects = []
            h = copy.deepcopy(mol.covalent_graph)
            h.remove_edge(*edge) #"breaks a bond" in the molecule graph
            connected_components = nx.algorithms.components.connected_components(h) #generates a set of nodes for each "fragment"
            for c in connected_components:

                subgraph = h.subgraph(c) #generates a subgraph from one set of nodes (this is a fragment graph)

                fragment_hash = weisfeiler_lehman_graph_hash( #saves the hash of this graph
                    subgraph, node_attr="specie"
                )

                fragment_hashes.append(fragment_hash) #adds each fragment hash to the fragment hash list

                if self.neighborhood_width is not None:
                    neighborhood_hashes = {}
                    for i in c:
                        hash_list = []
                        for d in range(self.neighborhood_width):
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
                    fragment_object = FragmentObject(
                        fragment_hash,
                        c,
                        neighborhood_hashes,
                        subgraph,
                        [i for i in c if i in edge[0:2]]
                    )
                    fragment_objects.append(fragment_object)

            fragment_complex = FragmentComplex(                                          #saves a FragmentComplex object after both fragment_hashes have been
                len(fragment_hashes), 1, [edge[0:2]], fragment_hashes, fragment_objects  #added to the list of fragments with len(fragments) fragments, 1 bond broken, the identity
            )                                                                            #of the bond broken (as a list containing one tuple), and the list of fragment hashes

            if len(fragment_hashes) == 1 and not self.allow_ring_opening:
                pass
            else:
                mol.fragment_data.append(fragment_complex) #append the above FragmentComplex object to the molecule's fragment_data list

        return False


class has_covalent_ring(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        # if mol is a metal, mol.covalent_graph is empty
        if mol.formula in m_formulas:
            mol.has_covalent_ring = False
        else:
            mol.has_covalent_ring = not nx.is_tree(mol.covalent_graph)

        if mol.has_covalent_ring:
            mol.ring_fragment_data = []

        return mol.has_covalent_ring


class covalent_ring_fragments(MSONable):
    # NOTE: haven't added fragment neighborhood hashes here yet!
    def __init__(self):
        pass

    def __call__(self, mol):
        # maps edge to graph with that edge removed
        ring_edges = {}

        for edge in mol.covalent_graph.edges:
            h = copy.deepcopy(mol.covalent_graph)
            h.remove_edge(*edge)
            if nx.is_connected(h):
                ring_edges[edge] = {
                    "modified_graph": h,
                    "node_set": set([edge[0], edge[1]]),
                }

        for ring_edge_1, ring_edge_2 in combinations(ring_edges, 2):

            if ring_edges[ring_edge_1]["node_set"].isdisjoint(
                ring_edges[ring_edge_2]["node_set"]
            ):

                potential_edges = [
                    (ring_edge_1[0], ring_edge_2[0], 0),
                    (ring_edge_1[0], ring_edge_2[1], 0),
                    (ring_edge_1[1], ring_edge_2[0], 0),
                    (ring_edge_1[1], ring_edge_2[1], 0),
                ]

                one_bond_away = False
                for ring_edge_3 in ring_edges:
                    if ring_edge_3 in potential_edges:
                        one_bond_away = True

                if one_bond_away:
                    h = copy.deepcopy(ring_edges[ring_edge_1]["modified_graph"])
                    h.remove_edge(*ring_edge_2)
                    if nx.is_connected(h):
                        continue
                    else:
                        fragments = []
                        connected_components = (
                            nx.algorithms.components.connected_components(h)
                        )
                        for c in connected_components:

                            subgraph = h.subgraph(c)

                            fragment_hash = weisfeiler_lehman_graph_hash(
                                subgraph, node_attr="specie"
                            )

                            fragments.append(fragment_hash)

                        fragment_complex = FragmentComplex(
                            len(fragments),
                            2,
                            [ring_edge_1[0:2], ring_edge_2[0:2]],
                            fragments,
                        )

                        mol.ring_fragment_data.append(fragment_complex)

        return False


class metal_complex(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        # if mol is a metal, it isn't a metal complex
        if mol.formula in m_formulas:
            return False

        return not nx.is_connected(mol.covalent_graph)


class h_atom_filter(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        # if mol is H+, H0, or H-
        return mol.formula == "H1"


class oh_plus_filter(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        # if mol is OH+
        return mol.formula == "H1 O1" and mol.charge == 1


class fix_hydrogen_bonding(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        if mol.num_atoms > 1:
            for i in range(mol.num_atoms):
                if mol.species[i] == "H":

                    adjacent_atoms = []

                    for bond in mol.graph.edges:
                        if i in bond[0:2]:

                            if i == bond[0]:
                                adjacent_atom = bond[1]
                            else:
                                adjacent_atom = bond[0]

                            displacement = (
                                mol.atom_locations[adjacent_atom]
                                - mol.atom_locations[i]
                            )

                            dist = np.inner(displacement, displacement)

                            adjacent_atoms.append((adjacent_atom, dist))

                    closest_atom, _ = min(adjacent_atoms, key=lambda pair: pair[1])

                    for adjacent_atom, _ in adjacent_atoms:
                        if adjacent_atom != closest_atom:
                            mol.graph.remove_edge(i, adjacent_atom)
                            if adjacent_atom in mol.covalent_graph:
                                mol.covalent_graph.remove_edge(i, adjacent_atom)

        return False


class bad_metal_coordination(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):

        if mol.formula not in m_formulas:

            if (
                len(metals.intersection(set(mol.species))) > 0
                and mol.number_of_coordination_bonds == 0
            ):

                return True

        return False


class set_solvation_correction(MSONable):
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

    def __init__(self, solvation_env):
        self.solvation_env = solvation_env

    def __call__(self, mol):
        correction = 0.0
        mol.number_of_coordination_bonds = 0

        for i in mol.m_inds:

            species = mol.species[i]
            partial_charge = mol.partial_charges_nbo[i]

            if partial_charge < 1.2:
                effective_charge = "_1"
            elif partial_charge >= 1.2:
                effective_charge = "_2"

            coordination_partners = list()
            species_charge = species + effective_charge
            radius = self.solvation_env["coordination_radius"][species_charge]

            for j in range(mol.num_atoms):
                if j != i:
                    displacement_vector = mol.atom_locations[j] - mol.atom_locations[i]
                    if np.inner(
                        displacement_vector, displacement_vector
                    ) < radius**2 and (
                        mol.partial_charges_resp[j] < 0
                        or mol.partial_charges_mulliken[j] < 0
                        or mol.partial_charges_nbo[j] < 0
                    ):
                        if not mol.graph.has_edge(i, j):
                            mol.graph.add_edge(i, j)
                        coordination_partners.append(j)

            number_of_coordination_bonds = len(coordination_partners)
            mol.number_of_coordination_bonds += number_of_coordination_bonds
            correction += self.solvation_env["solvation_correction"][species_charge] * (
                self.solvation_env["max_number_of_coordination_bonds"][species_charge]
                - number_of_coordination_bonds
            )

        mol.solvation_correction = correction
        return False


class species_default_true(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        return True


def compute_graph_hashes(mol):
    mol.total_hash = weisfeiler_lehman_graph_hash(mol.graph, node_attr="specie")

    mol.covalent_hash = weisfeiler_lehman_graph_hash(
        mol.covalent_graph, node_attr="specie"
    )

    return False


class neutral_metal_filter(MSONable):
    def __init__(self, cutoff):
        self.cutoff = cutoff

    def __call__(self, mol):

        for i in mol.m_inds:
            if mol.species[i] in metals and mol.partial_charges_nbo[i] < self.cutoff:
                return True

        return False


class charge_too_big(MSONable):
    def __init__(self):
        pass

    def __call__(self, mol):
        if mol.charge > 1 or mol.charge < -1:
            return True

        else:
            return False


# any species filter which modifies bonding has to come before
# any filter checking for connectivity (which includes the metal-centric complex filter)

li_species_decision_tree = [
    (fix_hydrogen_bonding(), Terminal.KEEP),
    (set_solvation_correction(li_ec), Terminal.KEEP),
    (charge_too_big(), Terminal.DISCARD),
    (neutral_metal_filter(0.1), Terminal.DISCARD),
    (compute_graph_hashes, Terminal.KEEP),
    (metal_ion_filter(), Terminal.DISCARD),
    (bad_metal_coordination(), Terminal.DISCARD),
    (mol_not_connected(), Terminal.DISCARD),
    (metal_complex(), Terminal.DISCARD),
    (spin_multiplicity_filter(0.4), Terminal.DISCARD),
    (add_star_hashes(), Terminal.KEEP),
    (add_unbroken_fragment(), Terminal.KEEP),
    (add_single_bond_fragments(), Terminal.KEEP),
    # (has_covalent_ring(), [
    #     (covalent_ring_fragments(), Terminal.KEEP),
    #     (species_default_true(), Terminal.KEEP)
    # ]),
    (species_default_true(), Terminal.KEEP),
]

mg_species_decision_tree = [
    (fix_hydrogen_bonding(), Terminal.KEEP),
    (set_solvation_correction(mg_g2), Terminal.KEEP),
    (neutral_metal_filter(0.5), Terminal.DISCARD),
    (compute_graph_hashes, Terminal.KEEP),
    (metal_ion_filter(), Terminal.DISCARD),
    (bad_metal_coordination(), Terminal.DISCARD),
    (mol_not_connected(), Terminal.DISCARD),
    (metal_complex(), Terminal.DISCARD),
    (add_star_hashes(), Terminal.KEEP),
    (add_unbroken_fragment(), Terminal.KEEP),
    (add_single_bond_fragments(), Terminal.KEEP),
    (species_default_true(), Terminal.KEEP),
]

nonmetal_species_decision_tree = [
    (fix_hydrogen_bonding(), Terminal.KEEP),
    (compute_graph_hashes, Terminal.KEEP),
    (add_star_hashes(), Terminal.KEEP),
    (add_unbroken_fragment(), Terminal.KEEP),
    (add_single_bond_fragments(allow_ring_opening=False), Terminal.KEEP),
    (species_default_true(), Terminal.KEEP),
]


euvl_species_decision_tree = [
    (fix_hydrogen_bonding(), Terminal.KEEP),
    (h_atom_filter(), Terminal.DISCARD),
    (oh_plus_filter(), Terminal.DISCARD),
    (compute_graph_hashes, Terminal.KEEP),
    (add_star_hashes(), Terminal.KEEP),
    (add_unbroken_fragment(neighborhood_width=3), Terminal.KEEP),
    (add_single_bond_fragments(allow_ring_opening=False, neighborhood_width=3), Terminal.KEEP),
    (species_default_true(), Terminal.KEEP),
]
