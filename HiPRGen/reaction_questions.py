import math
import numpy as np
from HiPRGen.mol_entry import MoleculeEntry
from enum import Enum
from functools import partial
from HiPRGen.constants import *
from scipy.optimize import linear_sum_assignment

"""
The reaction decision tree:

A question is a function q(reaction, mol_entries, params) -> Bool

reaction is a dict:

        reaction = {
            'reactants' : reactant indices
            'products' : product indices,
            'number_of_reactants',
            'number_of_products'}

params is a dict:


        params = {
           'temperature'
           'electron_free_energy'
        }

The lists of reactant and product indices always have length two. We use -1 when there is a only a single reactant or product.

The questions can also set reaction['rate'] and reaction['dG']

Questions will be writable by hand, or we could have machine learning filters.

A node is either a Terminal or a non empty list [(question, node)]

class Terminal(Enum):
    KEEP = 1
    DISCARD = -1

For the return value of a question, True means travel to this node and False means try next question in the list.

for non terminal nodes, it is an error if every question returns False. i.e getting stuck at a non terminal node is an error.

Once a Terminal node is reached, it tells us whether to keep or discard the reaction.

logging decision tree: The dispatcher takes a second decision tree as an argument, the logging decision tree. Reactions which return Terminal.KEEP from the logging decision tree will be logged in the generation report, with location specified by the argument generation_report_path


reaction atom mappings: compute_atom_mapping sets reaction['atom_mapping'] which is a dictionary sending (reactant_num, atom_num) to (product_num, atom_num).

"""

class Terminal(Enum):
    KEEP = 1
    DISCARD = -1

def run_decision_tree(
        reaction,
        mol_entries,
        params,
        decision_tree,
        decision_pathway=None):
    node = decision_tree

    while type(node) == list:
        next_node = None
        for (question, new_node) in node:
            if question(reaction, mol_entries, params):

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
        raise Exception(
            """
            unexpected node type reached.
            this is usually caused because none of the questions in some node returned True.
            """)



def default_rate(dG, params):
    kT = KB * params['temperature']
    max_rate = kT / PLANCK

    if dG < 0:
        rate = max_rate
    else:
        rate = max_rate * math.exp(- dG / kT)

    return rate

def dG_above_threshold(threshold, reaction, mol_entries, params):
    dG = 0.0

    # positive dCharge means electrons are lost
    dCharge = 0.0

    for i in range(reaction['number_of_reactants']):
        reactant_index = reaction['reactants'][i]
        mol = mol_entries[reactant_index]
        dG -= mol.get_free_energy()
        dCharge -= mol.charge

    for j in range(reaction['number_of_products']):
        product_index = reaction['products'][j]
        mol = mol_entries[product_index]
        dG += mol.get_free_energy()
        dCharge += mol.charge

    dG += dCharge * params['electron_free_energy']


    if dG > threshold:
        return True
    else:
        reaction['dG'] = dG
        reaction['rate'] = default_rate(dG, params)
        return False

def is_redox_reaction(reaction, mol_entries, params):
    # positive dCharge means electrons are lost
    dCharge = 0.0

    for i in range(reaction['number_of_reactants']):
        reactant_index = reaction['reactants'][i]
        mol = mol_entries[reactant_index]
        dCharge -= mol.charge

    for j in range(reaction['number_of_products']):
        product_index = reaction['products'][j]
        mol = mol_entries[product_index]
        dCharge += mol.charge

    if dCharge == 0:
        return False
    else:
        return True


def too_many_reactants_or_products(reaction, mols, params):
    if (reaction['number_of_reactants'] != 1 or
        reaction['number_of_products'] != 1):
        return True
    else:
        return False

def dcharge_too_large(reaction, mol_entries, params):
    dCharge = 0.0

    for i in range(reaction['number_of_reactants']):
        reactant_index = reaction['reactants'][i]
        mol = mol_entries[reactant_index]
        dCharge -= mol.charge

    for j in range(reaction['number_of_products']):
        product_index = reaction['products'][j]
        mol = mol_entries[product_index]
        dCharge += mol.charge

    if abs(dCharge) > 1:
        return True
    else:
        return False

def reactant_and_product_not_isomorphic(reaction, mols, params):
    reactant = mols[reaction['reactants'][0]]
    product = mols[reaction['products'][0]]
    if reactant.mol_graph.isomorphic_to(product.mol_graph):
        return False
    else:
        return True

def default_true(reaction, mols, params):
    return True


def is_A_B_to_A_C_where_A_not_metal_atom(reaction, mols, params):
    reactants_set = set([])
    products_set = set([])
    for i in range(reaction['number_of_reactants']):
        reactant_index = reaction['reactants'][i]
        reactants_set.add(reactant_index)

    for j in range(reaction['number_of_products']):
        product_index = reaction['products'][j]
        products_set.add(product_index)

    common_species = reactants_set.intersection(products_set)
    if len(common_species) == 0:
        return False
    else:
        mol = list(common_species)[0]
        if mols[mol].species == ['Li']:
            return False
        else:
            return True

def reaction_is_decomposable(reaction, mols, params):
    if reaction['number_of_reactants'] < 2 or reaction['number_of_products'] < 2:
        return False
    else:
        reactant_formulas = []
        product_formulas = []
        for i in range(reaction['number_of_reactants']):
            reactant = mols[reaction['reactants'][i]]
            reactant_formulas.append(reactant.formula)

        for i in range(reaction['number_of_products']):
            product = mols[reaction['products'][i]]
            product_formulas.append(product.formula)

        reactant_formulas_set = set(reactant_formulas)
        product_formulas_set = set(product_formulas)

        if len(reactant_formulas_set.intersection(product_formulas_set)) > 0:
            return True
        else:
            return False

def star_count_diff_above_threshold(threshold, reaction, mols, params):
    reactant_stars = {}
    product_stars = {}
    tags = set()

    for i in range(reaction['number_of_reactants']):
        reactant_index = reaction['reactants'][i]
        mol = mols[reactant_index]
        for h in mol.neighborhood_hashes[1].values():
            tags.add(h)
            if h in reactant_stars:
                reactant_stars[h] += 1
            else:
                reactant_stars[h] = 1

    for j in range(reaction['number_of_products']):
        product_index = reaction['products'][j]
        mol = mols[product_index]
        for h in mol.neighborhood_hashes[1].values():
            tags.add(h)
            if h in product_stars:
                product_stars[h] += 1
            else:
                product_stars[h] = 1

    count = 0

    for tag in tags:
        count += abs(reactant_stars.get(tag,0) - product_stars.get(tag,0))

    if count > threshold:
        return True
    else:
        return False



def compute_atom_mapping(radius_bound, reaction, mols, params):

    """
    atom mapping is in NP. It can be translated into a mixed integer programming problem, but this is too slow for HiPRGen. Instead, we use a combination of the heuristics described in "comparing heuristics gor graph edit distance computation, D.Blumenthal", to transform the atom mapping problem into a linear sum assignment problem that is solvable in linear time.

    the string mol.neighborhood_hashes[radius][atom_num] is an isomorphism preserving hash of the neighborhood arount atom_num with given radius. We want to map atoms in the reactants to atoms in the products with as much overlapping neighborhood as possible.
    """

    reactant_mapping = []
    product_mapping = []

    for i in range(reaction['number_of_reactants']):
        for j in range(mols[reaction['reactants'][i]].num_atoms):
            reactant_mapping.append((i,j))


    for i in range(reaction['number_of_products']):
        for j in range(mols[reaction['products'][i]].num_atoms):
            product_mapping.append((i,j))



    total_num_atoms = len(reactant_mapping)
    cost = np.zeros((total_num_atoms, total_num_atoms))

    for i in range(total_num_atoms):
        for j in range(total_num_atoms):
            reactant_num, reactant_atom_num = reactant_mapping[i]
            reactant_id = reaction['reactants'][reactant_num]
            reactant = mols[reactant_id]


            product_num, product_atom_num = product_mapping[j]
            product_id = reaction['products'][product_num]
            product = mols[product_id]

            bool_array = [
                reactant.neighborhood_hashes[r][reactant_atom_num] ==
                product.neighborhood_hashes[r][product_atom_num] for
                r in range(radius_bound) ]



            cost[i,j] = len([b for b in bool_array if b])


    row_ind, col_ind = linear_sum_assignment(cost, maximize=True)

    atom_mapping = dict(
        [ (reactant_mapping[row_ind[i]],
           (product_mapping[col_ind[i]], cost[row_ind[i], col_ind[i]]))
          for i in range(total_num_atoms) ])

    reaction['atom_mapping'] = atom_mapping

    return True


def atom_mapping_cost_above_threshold(threshold, reaction, mols, params):

    total_cost = 0.0
    atom_mapping = reaction['atom_mapping']
    for reactant_atom in atom_mapping:
        total_cost += atom_mapping[reactant_atom][1]

    return total_cost



standard_reaction_decision_tree = [
    (partial(dG_above_threshold, 0.5), Terminal.DISCARD),

    # redox branch
    (is_redox_reaction, [

        (too_many_reactants_or_products, Terminal.DISCARD),
        (dcharge_too_large, Terminal.DISCARD),

        # we may want to allow redox reactions between non
        # isomorphic species at some point
        (reactant_and_product_not_isomorphic, Terminal.DISCARD),
        (default_true, Terminal.KEEP)]),

    (reaction_is_decomposable, Terminal.DISCARD),

    # discard reactions of the form A+B->A+C unless A is a Li atom
    (is_A_B_to_A_C_where_A_not_metal_atom, Terminal.DISCARD),

    (partial(star_count_diff_above_threshold, 4), Terminal.DISCARD),

    (partial(compute_atom_mapping, 4),

     [
         (default_true, Terminal.KEEP)
     ])

    ]


# useful for testing file system performance
minimal_reaction_decision_tree = [
    (partial(dG_above_threshold, 0.5), Terminal.DISCARD),
    (default_true, Terminal.KEEP)
    ]

standard_logging_decision_tree = Terminal.DISCARD

def no_atom_mapping(reaction, mols, params):
    if 'atom_mapping' in reaction:
        return False
    else:
        return True

