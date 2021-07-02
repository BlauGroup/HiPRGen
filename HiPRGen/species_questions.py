from HiPRGen.mol_entry import MoleculeEntry
from enum import Enum
import networkx as nx
import copy



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

def add_stars(mol):
    species = mol.species
    bonds = mol.bonds

    for i in range(mol.num_atoms):
        center = species[i]
        boundary = {}

        for bond in bonds:
            if i in bond:
                if bond[0] == i:
                    end = species[bond[1]]
                if bond[1] == i:
                    end = species[bond[0]]

                if end in boundary:
                    boundary[end] += 1
                else:
                    boundary[end] = 1

        mol.stars[i] = (center, boundary)

    return False


def add_covalent_star_counts(mol):

    for i in mol.stars:
        center, boundary = mol.stars[i]

        if center not in metals:
            boundary_set = frozenset(
                [pair for pair in boundary.items() if pair[0] not in metals])

            tag = (center, boundary_set)

            if tag in mol.covalent_star_counts:
                mol.covalent_star_counts[tag] += 1
            else:
                mol.covalent_star_counts[tag] = 1

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
    # add_covalent_bond_count and add_covalent_stars always returns False
    # the associated Terminal nodes are never reached.
    (mol_not_connected, Terminal.DISCARD),
    (metal_ion_filter, Terminal.DISCARD),
    (metal_complex, Terminal.DISCARD),
    (add_stars, Terminal.KEEP),
    (add_covalent_star_counts, Terminal.KEEP),
    (default_true, Terminal.KEEP)
    ]
