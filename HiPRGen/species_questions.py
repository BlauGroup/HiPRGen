from HiPRGen.mol_entry import MoleculeEntry
from enum import Enum
import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
import copy
from functools import partial


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



metals = frozenset(["Li", "Na", "K", "Mg", "Ca", "Zn", "Al"])
m_formulas = frozenset([m + "1" for m in metals])

def metal_ion_filter(mol_entry):
    "only allow positively charged metal ions"
    if mol_entry.formula in m_formulas and mol_entry.charge <= 0:
        return True
    else:
        return False


def mol_not_connected(mol):
    return not nx.is_weakly_connected(mol.graph)


def add_star_hashes(mol):

    m_inds = [
        i for i, x in enumerate(mol.species) if x in metals
    ]

    g = copy.deepcopy(mol.graph.to_undirected())
    g.remove_nodes_from(m_inds)

    for i in range(mol.num_atoms):
        if i not in m_inds:
            neighborhood = nx.generators.ego.ego_graph(
                g,
                i,
                1,
                undirected=True)

            mol.star_hashes[i] = weisfeiler_lehman_graph_hash(
                neighborhood,
                node_attr='specie')

    return False

def add_total_hashes(mol):
    m_inds = [
        i for i, x in enumerate(mol.species) if x in metals
    ]

    g = copy.deepcopy(mol.graph.to_undirected())

    mol.total_hash = weisfeiler_lehman_graph_hash(
        g,
        node_attr='specie')

    g.remove_nodes_from(m_inds)

    mol.covalent_hash = weisfeiler_lehman_graph_hash(
        g,
        node_attr='specie')


    return False


def add_fragment_hashes(mol):
    if mol.formula in m_formulas:
        return False

    m_inds = [
        i for i, x in enumerate(mol.species) if x in metals
    ]

    g = copy.deepcopy(mol.graph.to_undirected())
    g.remove_nodes_from(m_inds)

    for edge in g.edges:
        h = copy.deepcopy(g)
        h.remove_edge(*edge)
        connected_components = nx.algorithms.components.connected_components(h)
        fragment_hashes = [
            weisfeiler_lehman_graph_hash(
                h.subgraph(c),
                node_attr='specie')

            for c in connected_components
            ]

        mol.fragment_hashes.append(fragment_hashes)

    return False



def metal_complex(mol):
    # if mol is a metal, it isn't a metal complex
    if mol.formula in m_formulas:
        return False

    # if mol has a metal, check if removing that metal disconnects.
    elif any([x in mol.formula for x in metals]):
        m_inds = [
            i for i, x in enumerate(mol.species) if x in metals
        ]
        g = copy.deepcopy(mol.mol_graph)
        g.remove_nodes(m_inds)
        return not nx.is_weakly_connected(g.graph)

    # no metal atoms
    else:
        return False

def default_true(mol):
    return True


standard_mol_decision_tree = [
    (mol_not_connected, Terminal.DISCARD),
    (metal_ion_filter, Terminal.DISCARD),
    (metal_complex, Terminal.DISCARD),
    (add_star_hashes, Terminal.KEEP),
    (add_total_hashes, Terminal.KEEP),
    (add_fragment_hashes, Terminal.KEEP),
    (default_true, Terminal.KEEP)
    ]
