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


def add_covalent_bond_count(mol_entry):
    mol_bond_count = {}
    species = mol_entry.species
    bonds = mol_entry.bonds
    for bond in bonds:
        species_0 = species[bond[0]]
        species_1 = species[bond[1]]
        tag = frozenset([species_0, species_1])
        if len(metals.intersection(tag)) == 0:
            if tag in mol_bond_count:
                mol_bond_count[tag] += 1
            else:
                mol_bond_count[tag] = 1

    mol_entry.aux_data['bond_count'] = mol_bond_count
    return False

def add_covalent_stars(mol_entry):
    species = mol_entry.species
    bonds = mol_entry.bonds
    stars = {}

    for i in range(mol_entry.num_atoms):
        if species[i] not in metals:
            end_counts = {}

            for bond in bonds:
                if i in bond:
                    if bond[0] == i:
                        end = species[bond[1]]
                    if bond[1] == i:
                        end = species[bond[0]]

                    if end in end_counts:
                        end_counts[end] += 1
                    else:
                        end_counts[end] = 1

            star = (species[i], frozenset(end_counts.items()))
            if star in stars:
                stars[star] += 1
            else:
                stars[star] = 1

    mol_entry.aux_data['stars'] = stars
    return False


def mol_not_connected(mol):
    return not nx.is_weakly_connected(mol.graph)


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
    (add_covalent_bond_count,Terminal.KEEP),
    (add_covalent_stars, Terminal.KEEP),
    (mol_not_connected, Terminal.DISCARD),
    (metal_complex, Terminal.DISCARD),
    (default_true, Terminal.KEEP)
    ]
