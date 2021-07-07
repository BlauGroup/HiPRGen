import math
from HiPRGen.mol_entry import MoleculeEntry
from enum import Enum
from functools import partial
from HiPRGen.constants import *

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
    if reactant.total_hash != product.total_hash:
        return True
    else:
        return False

def default_true(reaction, mols, params):
    return True


def star_count_diff_above_threshold(threshold, reaction, mols, params):
    reactant_stars = {}
    product_stars = {}
    tags = set()

    for i in range(reaction['number_of_reactants']):
        reactant_index = reaction['reactants'][i]
        mol = mols[reactant_index]
        for h in mol.star_hashes.values():
            tags.add(h)
            if h in reactant_stars:
                reactant_stars[h] += 1
            else:
                reactant_stars[h] = 1

    for j in range(reaction['number_of_products']):
        product_index = reaction['products'][j]
        mol = mols[product_index]
        for h in mol.star_hashes.values():
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

def reaction_is_covalent_decomposable(reaction, mols, params):
    if (reaction['number_of_reactants'] == 2 and
        reaction['number_of_products'] == 2):


        reactant_total_hashes = set()
        for i in range(reaction['number_of_reactants']):
            reactant_id = reaction['reactants'][i]
            reactant = mols[reactant_id]
            reactant_total_hashes.add(reactant.covalent_hash)

        product_total_hashes = set()
        for i in range(reaction['number_of_products']):
            product_id = reaction['products'][i]
            product = mols[product_id]
            product_total_hashes.add(product.covalent_hash)

        if len(reactant_total_hashes.intersection(product_total_hashes)) > 0:
            return True
        else:
            return False

    return False



def no_fragment_matching_found(reaction, mols, params):
    if reaction['number_of_reactants'] == 1:
        reactant = mols[reaction['reactants'][0]]
        reactant_fragments = set([frozenset([reactant.covalent_hash])])
        for fragments in reactant.fragment_hashes:
            reactant_fragments.add(
                frozenset(fragments))

    elif reaction['number_of_reactants'] == 2:
        reactant_0 = mols[reaction['reactants'][0]]
        reactant_1 = mols[reaction['reactants'][1]]

        reactant_fragments = set([
            frozenset([reactant_0.covalent_hash, reactant_1.covalent_hash])
        ])

        for fragments in reactant_0.fragment_hashes:
            reactant_fragments.add(
                frozenset(fragments + [reactant_1.covalent_hash]))

        for fragments in reactant_1.fragment_hashes:
            reactant_fragments.add(
                frozenset(fragments + [reactant_0.covalent_hash]))


    if reaction['number_of_products'] == 1:

        product = mols[reaction['products'][0]]
        product_fragments = set([frozenset([product.covalent_hash])])
        for fragments in product.fragment_hashes:
            product_fragments.add(
                frozenset(fragments))

    elif reaction['number_of_products'] == 2:

        product_0 = mols[reaction['products'][0]]
        product_1 = mols[reaction['products'][1]]

        product_fragments = set([
            frozenset([product_0.covalent_hash, product_1.covalent_hash])
        ])

        for fragments in product_0.fragment_hashes:
            product_fragments.add(
                frozenset(fragments + [product_1.covalent_hash]))

        for fragments in product_1.fragment_hashes:
            product_fragments.add(
                frozenset(fragments + [product_0.covalent_hash]))


    breakpoint()
    if len(reactant_fragments.intersection(product_fragments)) == 0:
        return True
    else:
        return False


standard_reaction_decision_tree = [
    (partial(dG_above_threshold, 0.5), Terminal.DISCARD),

    # redox branch
    (is_redox_reaction, [

        (too_many_reactants_or_products, Terminal.DISCARD),
        (dcharge_too_large, Terminal.DISCARD),
        (reactant_and_product_not_isomorphic, Terminal.DISCARD),
        (default_true, Terminal.KEEP)
    ]),

    (partial(star_count_diff_above_threshold, 4), Terminal.DISCARD),

    (reaction_is_covalent_decomposable, Terminal.DISCARD),

    (no_fragment_matching_found, Terminal.DISCARD),

    (default_true, Terminal.KEEP)
    ]


# useful for testing file system performance
minimal_reaction_decision_tree = [
    (partial(dG_above_threshold, 0.5), Terminal.DISCARD),
    (default_true, Terminal.KEEP)
    ]

standard_logging_decision_tree = Terminal.DISCARD

