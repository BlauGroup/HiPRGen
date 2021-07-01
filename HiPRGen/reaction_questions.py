import math
from HiPRGen.mol_entry import MoleculeEntry
from enum import Enum
from functools import partial
from HiPRGen.constants import *
from HiPRGen.atom_mapping import get_reaction_atom_mapping

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



def bond_count_diff_above_threshold(threshold, reaction, mol_entries, params):

    tags = set()

    for i in range(reaction['number_of_reactants']):
        reactant_index = reaction['reactants'][i]
        mol = mol_entries[reactant_index]
        tags.update(mol.aux_data['bond_count'].keys())

    for j in range(reaction['number_of_products']):
        product_index = reaction['products'][j]
        mol = mol_entries[product_index]
        tags.update(mol.aux_data['bond_count'].keys())

    count = 0

    for tag in tags:
        inter_count = 0
        for i in range(reaction['number_of_reactants']):
            reactant_index = reaction['reactants'][i]
            mol = mol_entries[reactant_index]
            inter_count += mol.aux_data['bond_count'].get(tag, 0)

        for j in range(reaction['number_of_products']):
            product_index = reaction['products'][j]
            mol = mol_entries[product_index]
            inter_count -= mol.aux_data['bond_count'].get(tag, 0)

        count += abs(inter_count)

    if count > threshold:
        return True
    else:
        return False

def star_count_diff_above_threshold(
        threshold,
        reaction,
        mol_entries,
        params):

    stars = set()

    for i in range(reaction['number_of_reactants']):
        reactant_index = reaction['reactants'][i]
        mol = mol_entries[reactant_index]
        stars.update(mol.aux_data['stars'].keys())

    for j in range(reaction['number_of_products']):
        product_index = reaction['products'][j]
        mol = mol_entries[product_index]
        stars.update(mol.aux_data['stars'].keys())


    count = 0

    for star in stars:
        inter_count = 0
        for i in range(reaction['number_of_reactants']):
            reactant_index = reaction['reactants'][i]
            mol = mol_entries[reactant_index]
            inter_count += mol.aux_data['stars'].get(star, 0)

        for j in range(reaction['number_of_products']):
            product_index = reaction['products'][j]
            mol = mol_entries[product_index]
            inter_count -= mol.aux_data['stars'].get(star, 0)

        count += abs(inter_count)

    if count > threshold:
        return True
    else:
        return False


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


def atom_mapping_branch(reaction, mols, params):
    reactant_mols_list = [mols[i] for i in reaction['reactants'] if i != -1]
    product_mols_list = [mols[i] for i in reaction['products'] if i != -1]

    l, r, i = get_reaction_atom_mapping(
        reactant_mols_list,
        product_mols_list)
    breakpoint()


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

    (partial(bond_count_diff_above_threshold, 2), Terminal.DISCARD),

    (partial(star_count_diff_above_threshold, 4), Terminal.DISCARD),

    # discard reactions of the form A+B->A+C unless A is a Li atom
    (is_A_B_to_A_C_where_A_not_metal_atom, Terminal.DISCARD),

    (default_true, Terminal.KEEP)
    ]


# useful for testing file system performance
minimal_reaction_decision_tree = [
    (partial(dG_above_threshold, 0.5), Terminal.DISCARD),
    (default_true, Terminal.KEEP)
    ]

standard_logging_decision_tree = Terminal.DISCARD
