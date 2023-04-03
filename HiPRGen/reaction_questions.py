import math
from HiPRGen.mol_entry import MoleculeEntry
from functools import partial
import itertools
import copy
import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from HiPRGen.constants import Terminal, ROOM_TEMP, KB, PLANCK, m_formulas
from monty.json import MSONable

"""
The reaction decision tree:

A question is a function q(reaction, mol_entries, params) -> Bool

reaction is a dict:

        reaction = { 'reactants' : reactant indices,
                     'products' : product indices,
                     'number_of_reactants',
                     'number_of_products'}
params is a dict:


        params = { 'temperature',
                   'electron_free_energy' }

The lists of reactant and product indices always have length two. We
use -1 when there is a only a single reactant or product.

The questions can also set reaction['rate'] and reaction['dG']

Questions will be writable by hand, or we could have machine learning
filters.

A node is either a Terminal or a non empty list [(question, node)]

class Terminal(Enum): KEEP = 1 DISCARD = -1

For the return value of a question, True means travel to this node and
False means try next question in the list.

for non terminal nodes, it is an error if every question returns
False. i.e getting stuck at a non terminal node is an error.

Once a Terminal node is reached, it tells us whether to keep or
discard the reaction.

logging decision tree: The dispatcher takes a second decision tree as
an argument, the logging decision tree. Reactions which return
Terminal.KEEP from the logging decision tree will be logged in the
generation report, with location specified by the argument
generation_report_path

"""

hydrogen_graph = nx.MultiGraph()
hydrogen_graph.add_node(0, specie="H")
hydrogen_hash = weisfeiler_lehman_graph_hash(hydrogen_graph, node_attr="specie")

fluorine_graph = nx.MultiGraph()
fluorine_graph.add_node(0, specie="F")
fluorine_hash = weisfeiler_lehman_graph_hash(fluorine_graph, node_attr="specie")

carbon_graph = nx.MultiGraph()
carbon_graph.add_node(0, specie="C")
carbon_hash = weisfeiler_lehman_graph_hash(carbon_graph, node_attr="specie")


def run_decision_tree(
    reaction, mol_entries, params, decision_tree, decision_pathway=None
):
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


def default_rate(dG_barrier, params):
    kT = KB * params["temperature"]
    max_rate = kT / PLANCK
    rate = max_rate * math.exp(-dG_barrier / kT)
    return rate


class dG_above_threshold(MSONable):
    def __init__(self, threshold, free_energy_type, constant_barrier, barrier_factor=0):

        self.threshold = threshold
        self.free_energy_type = free_energy_type
        self.constant_barrier = constant_barrier
        self.barrier_factor = barrier_factor

        if free_energy_type == "free_energy":
            self.get_free_energy = lambda mol, temperature: mol.get_free_energy(temperature)
        elif free_energy_type == "solvation_free_energy":
            self.get_free_energy = lambda mol, temperature: mol.solvation_correction + mol.get_free_energy(temperature)
        else:
            raise Exception("unrecognized free energy type")

    def __str__(self):
        return self.free_energy_type + " dG is above threshold=" + str(self.threshold)

    def __call__(self, reaction, mol_entries, params):

        dG = 0.0

        # positive dCharge means electrons are lost
        dCharge = 0.0

        for i in range(reaction["number_of_reactants"]):
            reactant_index = reaction["reactants"][i]
            mol = mol_entries[reactant_index]
            dG -= self.get_free_energy(mol, params["temperature"])
            dCharge -= mol.charge

        for j in range(reaction["number_of_products"]):
            product_index = reaction["products"][j]
            mol = mol_entries[product_index]
            dG += self.get_free_energy(mol, params["temperature"])
            dCharge += mol.charge

        dG += dCharge * params["electron_free_energy"]

        if dG > self.threshold:
            reaction["dG"] = dG
            if self.barrier_factor == 0:
                barrier = self.constant_barrier
            else:
                barrier = reaction["dG"] * self.barrier_factor
            reaction["dG_barrier"] = barrier
            reaction["rate"] = default_rate(barrier, params)
            return True
        else:
            reaction["dG"] = dG
            if dG < 0:
                barrier = self.constant_barrier
            else:
                barrier = dG + self.constant_barrier

            reaction["dG_barrier"] = barrier
            reaction["rate"] = default_rate(barrier, params)
            return False


class dG_below_threshold(MSONable):
    def __init__(self, threshold, free_energy_type, constant_barrier):

        self.threshold = threshold
        self.free_energy_type = free_energy_type
        self.constant_barrier = constant_barrier

        if free_energy_type == "free_energy":
            self.get_free_energy = lambda mol, temperature: mol.get_free_energy(temperature)
        elif free_energy_type == "solvation_free_energy":
            self.get_free_energy = lambda mol, temperature: mol.solvation_correction + mol.get_free_energy(temperature)
        else:
            raise Exception("unrecognized free energy type")

    def __str__(self):
        return self.free_energy_type + " dG is below threshold=" + str(self.threshold)

    def __call__(self, reaction, mol_entries, params):

        dG = 0.0

        # positive dCharge means electrons are lost
        dCharge = 0.0

        for i in range(reaction["number_of_reactants"]):
            reactant_index = reaction["reactants"][i]
            mol = mol_entries[reactant_index]
            dG -= self.get_free_energy(mol, params["temperature"])
            dCharge -= mol.charge

        for j in range(reaction["number_of_products"]):
            product_index = reaction["products"][j]
            mol = mol_entries[product_index]
            dG += self.get_free_energy(mol, params["temperature"])
            dCharge += mol.charge

        dG += dCharge * params["electron_free_energy"]

        if dG < self.threshold:
            return True
        else:
            return False


class is_redox_reaction(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "is redox reaction"

    def __call__(self, reaction, mol_entries, params):
        # positive dCharge means electrons are lost
        dCharge = 0.0

        for i in range(reaction["number_of_reactants"]):
            reactant_index = reaction["reactants"][i]
            mol = mol_entries[reactant_index]
            dCharge -= mol.charge

        for j in range(reaction["number_of_products"]):
            product_index = reaction["products"][j]
            mol = mol_entries[product_index]
            dCharge += mol.charge

        if dCharge == 0:
            reaction["is_redox"] = False
            return False
        else:
            reaction["is_redox"] = True
            return True


class add_electron_species(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "add electron species"

    def __call__(self, reaction, mol_entries, params):
        if "electron_species" in params:
            # positive dCharge means electrons are lost
            dCharge = 0.0

            if reaction["number_of_reactants"] != 1:
                return False

            if reaction["number_of_products"] != 1:
                return False

            dCharge -= mol_entries[reaction["reactants"][0]].charge
            dCharge += mol_entries[reaction["products"][0]].charge

            if dCharge > 0:
                reaction["products"] = (
                    reaction["products"][0],
                    params["electron_species"],
                )
                reaction["number_of_products"] = 2
            elif dCharge < 0:
                reaction["reactants"] = (
                    reaction["reactants"][0],
                    params["electron_species"],
                )
                reaction["number_of_reactants"] = 2
            else:
                pass
            return False
        else:
            return False


class too_many_reactants_or_products(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "too many reactants or products"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] != 1 or reaction["number_of_products"] != 1:
            return True
        else:
            return False


class more_than_one_reactant(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "more than one reactant"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] != 1:
            return True
        else:
            return False


class only_one_product(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "only one product"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_products"] == 1:
            return True
        else:
            return False


class metal_metal_reaction(MSONable):
    def __init__(self):
        pass

    def __call__(self, reaction, mol_entries, params):
        if (
            reaction["number_of_reactants"] == 1
            and reaction["number_of_products"] == 1
            and (mol_entries[reaction["reactants"][0]].formula in m_formulas)
            and (mol_entries[reaction["products"][0]].formula in m_formulas)
        ):

            return True
        else:
            return False


class dcharge_too_large(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "change in charge is too large"

    def __call__(self, reaction, mol_entries, params):
        dCharge = 0.0

        for i in range(reaction["number_of_reactants"]):
            reactant_index = reaction["reactants"][i]
            mol = mol_entries[reactant_index]
            dCharge -= mol.charge

        for j in range(reaction["number_of_products"]):
            product_index = reaction["products"][j]
            mol = mol_entries[product_index]
            dCharge += mol.charge

        if abs(dCharge) > 1:
            return True
        else:
            return False


def marcus_barrier(reaction, mol_entries, params):

    """
        Okay, so Marcus Theory.The math works out like so.∆G* = λ/4 (1 +
    ∆G / λ)^2 ∆G is the Gibbs free energy of the reaction, ∆G* is the
    energy barrier, and λ is the “reorganization energy” (basically the
    energy penalty for reorganizing the solvent environment to accommodate
    the change in local charge).The reorganization energy can be broken up
    into two terms, an inner term (“i”) representing the contribution from
    the first solvation shell and an outer term (“o”) representing the
    contribution from the bulk solvent: λ = λi + λoλo = ∆e/(8 pi ε0) (1/r
    - 1/R) (1/n^2 - 1/ε) where ∆e is the change in charge in terms of
    fundamental charge (1.602 * 10 ^-19 C), ε0 is the vacuum permittivity
    (8.854 * 10 ^-12 F/m), r is the first solvation shell radius (I
    usually just pick a constant, say 6 Angstrom), R is the distance to
    the electrode (again, for these purposes, just pick something - say
    7.5 Angstrom), n is the index of refraction (1.415 for EC) and ε is
    the relative dielectric (18.5 for EC/EMC).
    """

    reactant = mol_entries[reaction["reactants"][0]]
    product = mol_entries[reaction["products"][0]]
    dCharge = product.charge - reactant.charge
    n = 1.415  # index of refraction; variable
    eps = 18.5  # dielectric constant; variable

    r = 6.0  # in Angstrom
    R = 7.5  # in Angstrom

    eps_0 = 8.85419 * 10**-12  # vacuum permittivity
    e = 1.602 * 10**-19  # fundamental charge

    l_outer = e / (8 * math.pi * eps_0)
    l_outer *= (
        1 / r - 1 / (2 * R)
    ) * 10**10  # Converting to SI units; factor of 2 is because of different definitions of the distance to electrode
    l_outer *= 1 / n**2 - 1 / eps

    if dCharge == -1:
        vals = [reactant.electron_affinity, product.ionization_energy]
        vals_filtered = [v for v in vals if v is not None]
        l_inner = sum(vals_filtered) / len(vals_filtered)

    if dCharge == 1:
        vals = [reactant.ionization_energy, product.electron_affinity]
        vals_filtered = [v for v in vals if v is not None]
        l_inner = sum(vals_filtered) / len(vals_filtered)

    if l_inner < 0:
        l_inner = 0

    l = l_inner + l_outer

    dG = (
        product.free_energy
        - reactant.free_energy
        + dCharge * params["electron_free_energy"]
    )
    dG_barrier = l / 4 * (1 + dG / l) ** 2
    reaction["marcus_barrier"] = dG_barrier
    return False


class reactant_and_product_not_isomorphic(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "reactants and products are not covalent isomorphic"

    def __call__(self, reaction, mol_entries, params):
        reactant = mol_entries[reaction["reactants"][0]]
        product = mol_entries[reaction["products"][0]]
        if reactant.covalent_hash != product.covalent_hash:
            return True
        else:
            return False


class reaction_default_true(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "default true"

    def __call__(self, reaction, mol_entries, params):
        return True


class star_count_diff_above_threshold(MSONable):
    """
    if you want to filter out break-one-form-one reactions, the
    correct value for the threshold is 6.
    """

    def __init__(self, threshold):
        self.threshold = threshold

    def __str__(self):
        return "star count diff above threshold=" + str(self.threshold)

    def __call__(self, reaction, mol_entries, params):
        reactant_stars = {}
        product_stars = {}
        tags = set()

        for i in range(reaction["number_of_reactants"]):
            reactant_index = reaction["reactants"][i]
            mol = mol_entries[reactant_index]
            for h in mol.star_hashes.values():
                tags.add(h)
                if h in reactant_stars:
                    reactant_stars[h] += 1
                else:
                    reactant_stars[h] = 1

        for j in range(reaction["number_of_products"]):
            product_index = reaction["products"][j]
            mol = mol_entries[product_index]
            for h in mol.star_hashes.values():
                tags.add(h)
                if h in product_stars:
                    product_stars[h] += 1
                else:
                    product_stars[h] = 1

        count = 0

        for tag in tags:
            count += abs(reactant_stars.get(tag, 0) - product_stars.get(tag, 0))

        if count > self.threshold:
            return True
        else:
            return False


class reaction_is_charge_transfer(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "reaction is charge transfer"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] == 2 and reaction["number_of_products"] == 2:

            reactant_total_hashes = set()
            for i in range(reaction["number_of_reactants"]):
                reactant_id = reaction["reactants"][i]
                reactant = mol_entries[reactant_id]
                reactant_total_hashes.add(reactant.covalent_hash)

            product_total_hashes = set()
            for i in range(reaction["number_of_products"]):
                product_id = reaction["products"][i]
                product = mol_entries[product_id]
                product_total_hashes.add(product.covalent_hash)

            if len(reactant_total_hashes.intersection(product_total_hashes)) == 2:
                return True
            else:
                return False

        return False


class reaction_is_covalent_decomposable(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "reaction is covalent decomposable"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] == 2 and reaction["number_of_products"] == 2:

            reactant_total_hashes = set()
            for i in range(reaction["number_of_reactants"]):
                reactant_id = reaction["reactants"][i]
                reactant = mol_entries[reactant_id]
                reactant_total_hashes.add(reactant.covalent_hash)

            product_total_hashes = set()
            for i in range(reaction["number_of_products"]):
                product_id = reaction["products"][i]
                product = mol_entries[product_id]
                product_total_hashes.add(product.covalent_hash)

            if len(reactant_total_hashes.intersection(product_total_hashes)) > 0:
                return True
            else:
                return False

        return False


class reaction_is_radical_separation(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "reaction is radical separation"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] == 1 and reaction["number_of_products"] == 2:
            reactant_spin = mol_entries[reaction["reactants"][0]].spin_multiplicity
            if reactant_spin == 1:
                prod0_spin = mol_entries[reaction["products"][0]].spin_multiplicity
                if prod0_spin != 1:
                    return True
        return False


class reaction_is_charge_separation(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "reaction is charge separation"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] == 1 and reaction["number_of_products"] == 2:
            reactant_charge = mol_entries[reaction["reactants"][0]].charge
            prod0_charge = mol_entries[reaction["products"][0]].charge
            prod1_charge = mol_entries[reaction["products"][1]].charge
            if abs(prod0_charge) > abs(reactant_charge) or abs(prod1_charge) > abs(reactant_charge):
                return True
        elif reaction["number_of_reactants"] == 2 and reaction["number_of_products"] == 2:
            reactant0_charge = mol_entries[reaction["reactants"][0]].charge
            reactant1_charge = mol_entries[reaction["reactants"][1]].charge
            prod0_charge = mol_entries[reaction["products"][0]].charge
            prod1_charge = mol_entries[reaction["products"][1]].charge
            if reactant0_charge == 0 and reactant1_charge == 0 and (abs(prod0_charge) > 0 or abs(prod1_charge)) > 0:
                return True
        return False


class reactants_are_both_anions_or_both_cations(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "reactants are both anions or both cations"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] == 2:
            reactant0_charge = mol_entries[reaction["reactants"][0]].charge
            reactant1_charge = mol_entries[reaction["reactants"][1]].charge
            if reactant0_charge > 0 and reactant1_charge > 0:
                return True
            elif reactant0_charge < 0 and reactant1_charge < 0:
                return True
        return False


class two_closed_shell_reactants_and_two_open_shell_products(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "two closed shell reactants and two open shell products"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] == 2 and reaction["number_of_products"] == 2:
            reactant0_spin = mol_entries[reaction["reactants"][0]].spin_multiplicity
            reactant1_spin = mol_entries[reaction["reactants"][1]].spin_multiplicity
            product0_spin = mol_entries[reaction["products"][0]].spin_multiplicity
            product1_spin = mol_entries[reaction["products"][1]].spin_multiplicity
            if reactant0_spin == 1 and reactant1_spin == 1 and product0_spin > 1 and product1_spin > 1:
                return True
        return False


class metal_coordination_passthrough(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "metal coordination passthrough"

    def __call__(self, reaction, mol_entries, params):

        for i in range(reaction["number_of_reactants"]):
            reactant_id = reaction["reactants"][i]
            reactant = mol_entries[reactant_id]
            if reactant.formula in m_formulas:
                return True

        for i in range(reaction["number_of_products"]):
            product_id = reaction["products"][i]
            product = mol_entries[product_id]
            if product.formula in m_formulas:
                return True

        return False


class compositions_preclude_h_transfer(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "compositions preclude h transfer"

    def __call__(self, reaction, mol_entries, params):
        reactant_compositions = []
        reactant_charges = []
        for i in range(reaction["number_of_reactants"]):
            reactant_id = reaction["reactants"][i]
            reactant = mol_entries[reactant_id]
            reactant_compositions.append(reactant.molecule.composition)
            reactant_charges.append(reactant.molecule.charge)
            
        product_compositions = []
        product_charges = []
        for i in range(reaction["number_of_products"]):
            product_id = reaction["products"][i]
            product = mol_entries[product_id]
            product_compositions.append(product.molecule.composition)
            product_charges.append(product.molecule.charge)

        if len(reactant_compositions) != 2 or len(product_compositions) != 2:
            return True

        h_transfer_possible = None

        try:
            comp_diff = reactant_compositions[0] - product_compositions[0]
            if comp_diff.alphabetical_formula == "H1":
                if abs(reactant_charges[0] - product_charges[0]) > 1:
                    h_transfer_possible = False
                else:
                    h_transfer_possible = True
        except ValueError:
            try:
                comp_diff = reactant_compositions[1] - product_compositions[0]
                if comp_diff.alphabetical_formula == "H1":
                    if abs(reactant_charges[1] - product_charges[0]) > 1:
                        h_transfer_possible = False
                    else:
                        h_transfer_possible = True
            except ValueError:
                try:
                    comp_diff = reactant_compositions[1] - product_compositions[1]
                    if comp_diff.alphabetical_formula == "H1":
                        if abs(reactant_charges[1] - product_charges[1]) > 1:
                            h_transfer_possible = False
                        else:
                            h_transfer_possible = True
                except ValueError:
                    try:
                        comp_diff = reactant_compositions[0] - product_compositions[1]
                        if comp_diff.alphabetical_formula == "H1":
                            if abs(reactant_charges[0] - product_charges[1]) > 1:
                                h_transfer_possible = False
                            else:
                                h_transfer_possible = True
                    except ValueError:
                        h_transfer_possible = False

        return not h_transfer_possible


class fragment_matching_found(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "fragment matching found"

    def __call__(self, reaction, mol_entries, params):

        reactant_fragment_indices_list = []
        product_fragment_indices_list = []

        if reaction["number_of_reactants"] == 1: #creates a list of the indicies pointing to FragmentComplex objects
            reactant = mol_entries[reaction["reactants"][0]] #reactant is a mol_entry
            for i in range(len(reactant.fragment_data)):  #fragment_data is a list of FragmentComplex objects, where each
                reactant_fragment_indices_list.append([i]) #FragmentComplex object is basically a dictionary with four keys

        if reaction["number_of_reactants"] == 2: 
            reactant_0 = mol_entries[reaction["reactants"][0]]
            reactant_1 = mol_entries[reaction["reactants"][1]]
            for i in range(len(reactant_0.fragment_data)): #for each fragment of one reactant
                for j in range(len(reactant_1.fragment_data)): #look at each fragment of the other reactant
                    if (                                                    #true only when adding fragments of one reactant with the other 
                        reactant_0.fragment_data[i].number_of_bonds_broken  #unfragmented reactant?
                        + reactant_1.fragment_data[j].number_of_bonds_broken
                        <= 1 
                    ): 

                        reactant_fragment_indices_list.append([i, j]) #append a list to the list containing fragment indicies for both reactants

        if reaction["number_of_products"] == 1: #repeat for product indicies
            product = mol_entries[reaction["products"][0]]
            for i in range(len(product.fragment_data)):
                product_fragment_indices_list.append([i])

        if reaction["number_of_products"] == 2:
            product_0 = mol_entries[reaction["products"][0]]
            product_1 = mol_entries[reaction["products"][1]]
            for i in range(len(product_0.fragment_data)):
                for j in range(len(product_1.fragment_data)):
                    if (
                        product_0.fragment_data[i].number_of_bonds_broken
                        + product_1.fragment_data[j].number_of_bonds_broken
                        <= 1
                    ):

                        product_fragment_indices_list.append([i, j])

        viable_fragment_matches = []
        for reactant_fragment_indices in reactant_fragment_indices_list: #iterating over all reactant and product fragment indicies
            for product_fragment_indices in product_fragment_indices_list:
                reactant_fragment_count = 0
                product_fragment_count = 0
                reactant_bonds_broken = []
                product_bonds_broken = []

                reactant_hashes = dict()
                for reactant_index, frag_complex_index in enumerate(reactant_fragment_indices):
                    fragment_complex = mol_entries[                  #pulls out a fragment_complex whose index matches the above
                        reaction["reactants"][reactant_index]
                    ].fragment_data[frag_complex_index]

                    for bond in fragment_complex.bonds_broken:       #save what bonds are broken in this complex to reactant_bonds_broken
                        reactant_bonds_broken.append(
                            [(reactant_index, x) for x in bond] #first element of tuple is which reactant, x is a integer, bond is a tuple containing two numbers denoting the edge of a molecule graph
                        )

                    for i in range(fragment_complex.number_of_fragments):
                        reactant_fragment_count += 1
                        tag = fragment_complex.fragment_hashes[i]
                        if tag in reactant_hashes:
                            reactant_hashes[tag] += 1
                        else:
                            reactant_hashes[tag] = 1

                product_hashes = dict()
                for product_index, frag_complex_index in enumerate(
                    product_fragment_indices
                ):

                    fragment_complex = mol_entries[
                        reaction["products"][product_index]
                    ].fragment_data[frag_complex_index]

                    for bond in fragment_complex.bonds_broken:
                        product_bonds_broken.append([(product_index, x) for x in bond])

                    for i in range(fragment_complex.number_of_fragments):
                        product_fragment_count += 1
                        tag = fragment_complex.fragment_hashes[i]
                        if tag in product_hashes:
                            product_hashes[tag] += 1
                        else:
                            product_hashes[tag] = 1

                # don't consider fragmentations with both a ring opening and closing
                if (
                    reaction["number_of_reactants"] == 2
                    and reaction["number_of_products"] == 2
                    and reactant_fragment_count == 2
                    and product_fragment_count == 2
                ):
                    continue

                if reactant_hashes == product_hashes:
                    if hydrogen_hash in reactant_hashes:
                        reaction["reactant_bonds_broken"] = reactant_bonds_broken
                        reaction["product_bonds_broken"] = product_bonds_broken
                        reaction["hashes"] = reactant_hashes
                        reaction["reactant_fragment_count"] = reactant_fragment_count
                        reaction["product_fragment_count"] = product_fragment_count
                        return True
                    else:
                        tmp = {}
                        tmp["reactant_bonds_broken"] = reactant_bonds_broken
                        tmp["product_bonds_broken"] = product_bonds_broken
                        tmp["hashes"] = reactant_hashes
                        tmp["reactant_fragment_count"] = reactant_fragment_count
                        tmp["product_fragment_count"] = product_fragment_count
                        viable_fragment_matches.append(tmp)

        if len(viable_fragment_matches) > 0:
            min_frag_size = 1000000000
            if len(viable_fragment_matches) == 1:
                best_matching = viable_fragment_matches[0]
            else:
                for viable_match in viable_fragment_matches:
                    if len(viable_match["reactant_bonds_broken"]) == 0:
                        for reactant_index in reaction["reactants"]:
                            reactant = mol_entries[reactant_index]
                            if len(reactant.molecule) < min_frag_size:
                                min_frag_size = len(reactant.molecule)
                                best_matching = copy.deepcopy(viable_match)
                    else:
                        for l in viable_match["reactant_bonds_broken"]:
                            hot_reactant = mol_entries[reaction["reactants"][l[0][0]]]
                            hot_reactant_graph = copy.deepcopy(hot_reactant.covalent_graph)
                            edge = (l[0][1],l[1][1])
                            hot_reactant_graph.remove_edge(*edge)
                            connected_components = nx.algorithms.components.connected_components(hot_reactant_graph)
                            for c in connected_components:
                                subgraph = hot_reactant_graph.subgraph(c)
                                # if subgraph.number_of_nodes() == 1:
                                #     for node in subgraph:
                                #         print(nx.get_node_attributes(hot_reactant_graph, "specie")[node])
                                if subgraph.number_of_nodes() < min_frag_size:
                                    min_frag_size = subgraph.number_of_nodes()
                                    best_matching = copy.deepcopy(viable_match)
            reaction["reactant_bonds_broken"] = best_matching["reactant_bonds_broken"]
            reaction["product_bonds_broken"] = best_matching["product_bonds_broken"]
            reaction["hashes"] = best_matching["hashes"]
            reaction["reactant_fragment_count"] = best_matching["reactant_fragment_count"]
            reaction["product_fragment_count"] = best_matching["product_fragment_count"]
            return True

        return False


class single_reactant_single_product_not_atom_transfer(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "not intramolecular hydrogen/fluorine transfer"

    def __call__(self, reaction, mol_entries, params):
        if (
            reaction["number_of_reactants"] == 1
            and reaction["number_of_products"] == 1
            and len(reaction["reactant_bonds_broken"]) == 1
            and len(reaction["product_bonds_broken"]) == 1
            and hydrogen_hash not in reaction["hashes"]
            and fluorine_hash not in reaction["hashes"]
        ):
            return True

        return False


class not_h_transfer(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "not intermolecular hydrogen transfer"

    def __call__(self, reaction, mol_entries, params):
        if hydrogen_hash not in reaction["hashes"]:
            return True

        return False


class single_reactant_single_product(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "single reactant single product"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] == 1 and reaction["number_of_products"] == 1:
            return True

        return False


class single_reactant_double_product_ring_close(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "ring close"

    def __call__(self, reaction, mol_entries, params):

        if (
            reaction["number_of_reactants"] == 1
            and reaction["number_of_products"] == 2
            and len(reaction["reactant_bonds_broken"]) == 1
            and len(reaction["product_bonds_broken"]) == 1
            and reaction["product_fragment_count"] == 2
        ):

            return True

        return False


class h_abstraction_from_closed_shell_reactant(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "h abstraction from closed shell reactant"

    def __call__(self, reaction, mol_entries, params):

        if reaction["number_of_reactants"] == 2 and reaction["number_of_products"] == 2 and hydrogen_hash in reaction["hashes"]:
            hot_reactant_ind = reaction["reactant_bonds_broken"][0][0][0]
            cold_reactant_ind = 0
            if hot_reactant_ind == 0:
                cold_reactant_ind = 1
            hot_product_ind = reaction["product_bonds_broken"][0][0][0]
            hot_reactant = mol_entries[reaction["reactants"][hot_reactant_ind]]
            cold_reactant = mol_entries[reaction["reactants"][cold_reactant_ind]]
            hot_product = mol_entries[reaction["products"][hot_product_ind]]
            if hot_reactant.spin_multiplicity == 1:
                if hot_product.charge - cold_reactant.charge == 0:
                    return True

        return False


class h_minus_abstraction(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "h minus abstraction"

    def __call__(self, reaction, mol_entries, params):
        if reaction["number_of_reactants"] == 2 and reaction["number_of_products"] == 2 and hydrogen_hash in reaction["hashes"]:
            hot_reactant_ind = reaction["reactant_bonds_broken"][0][0][0]
            cold_reactant_ind = 0
            if hot_reactant_ind == 0:
                cold_reactant_ind = 1
            hot_product_ind = reaction["product_bonds_broken"][0][0][0]
            cold_reactant = mol_entries[reaction["reactants"][cold_reactant_ind]]
            hot_product = mol_entries[reaction["products"][hot_product_ind]]
            if hot_product.charge - cold_reactant.charge == -1:
                return True

        return False


class concerted_metal_coordination(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "concerted metal coordination"

    def __call__(self, reaction, mol_entries, params):

        if reaction["number_of_reactants"] == 2 and reaction["number_of_products"] == 2:

            reactant_0 = mol_entries[reaction["reactants"][0]]
            reactant_1 = mol_entries[reaction["reactants"][1]]
            product_0 = mol_entries[reaction["products"][0]]
            product_1 = mol_entries[reaction["products"][1]]

            if (
                reactant_0.formula in m_formulas
                or reactant_1.formula in m_formulas
                or product_0.formula in m_formulas
                or product_1.formula in m_formulas
            ):
                return True
            else:
                return False

        return False


class concerted_metal_coordination_one_product(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "concerted metal coordination one product"

    def __call__(self, reaction, mol_entries, params):

        if reaction["number_of_reactants"] == 2 and reaction["number_of_products"] == 1:

            reactant_0 = mol_entries[reaction["reactants"][0]]
            reactant_1 = mol_entries[reaction["reactants"][1]]
            product = mol_entries[reaction["products"][0]]

            reactant_covalent_hashes = set(
                [reactant_0.covalent_hash, reactant_1.covalent_hash]
            )

            if (
                reactant_0.formula in m_formulas or reactant_1.formula in m_formulas
            ) and product.covalent_hash not in reactant_covalent_hashes:
                return True
            else:
                return False

        return False


class concerted_metal_coordination_one_reactant(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "concerted metal coordination one reactant"

    def __call__(self, reaction, mol_entries, params):

        if reaction["number_of_reactants"] == 1 and reaction["number_of_products"] == 2:

            product_0 = mol_entries[reaction["products"][0]] #use this to get a mol entry
            product_1 = mol_entries[reaction["products"][1]]
            reactant = mol_entries[reaction["reactants"][0]]

            product_covalent_hashes = set(
                [product_0.covalent_hash, product_1.covalent_hash]
            )

            if (
                product_0.formula in m_formulas or product_1.formula in m_formulas
            ) and reactant.covalent_hash not in product_covalent_hashes:
                return True
            else:
                return False

        return False


class single_reactant_with_ring_break_two(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "single reactant with a ring, break two"

    def __call__(self, reaction, mol_entries, params):
        if (
            reaction["number_of_reactants"] == 1
            and reaction["number_of_products"] == 2
            and mol_entries[reaction["reactants"][0]].has_covalent_ring
        ):

            reactant = mol_entries[reaction["reactants"][0]]
            product_1 = mol_entries[reaction["products"][0]]
            product_2 = mol_entries[reaction["products"][1]]
            for fragment_complex in reactant.ring_fragment_data:
                if set(fragment_complex.fragment_hashes) == set(
                    [product_1.covalent_hash, product_2.covalent_hash]
                ):
                    return True

        return False


class single_product_with_ring_form_two(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "single product with a ring, form two"

    def __call__(self, reaction, mol_entries, params):
        if (
            reaction["number_of_reactants"] == 2
            and reaction["number_of_products"] == 1
            and mol_entries[reaction["products"][0]].has_covalent_ring
        ):

            product = mol_entries[reaction["products"][0]]
            reactant_1 = mol_entries[reaction["reactants"][0]]
            reactant_2 = mol_entries[reaction["reactants"][1]]
            for fragment_complex in product.ring_fragment_data:
                if set(fragment_complex.fragment_hashes) == set(
                    [reactant_1.covalent_hash, reactant_2.covalent_hash]
                ):
                    return True

        return False


class reaction_is_hindered(MSONable):
    def __init__(self):
        pass

    def __str__(self):
        return "reaction is hindered"

    def __call__(self, reaction, mol_entries, params):

        hot_reactant_atoms = []

        for l in reaction["reactant_bonds_broken"]: #finds the indicies for the atoms in the broken bond
            for t in l:
                hot_reactant = mol_entries[reaction["reactants"][t[0]]]
                hot_reactant_atoms.append(t[1])

        hot_product_atoms = []

        for l in reaction["product_bonds_broken"]: #finds the indicies for the atoms in the formed bond
            for t in l:
                hot_product = mol_entries[reaction["products"][t[0]]]
                hot_product_atoms.append(t[1])

        reaction_methyl_test = []
        reactant_num_carbon_neighbors = 0
        for atom in hot_reactant_atoms:
            reactant_num_hydrogens = 0
            if hot_reactant.mol_graph.get_coordination_of_site(atom) == 4: #only care about sp3 hybidized carbons
                neighbor_list = hot_reactant.mol_graph.get_connected_sites(atom)
                for neighbor in neighbor_list:
                    neighbor_index = neighbor[2]
                    if hot_reactant.mol_graph.get_coordination_of_site(neighbor_index) == 4: #if neighbor is also sp3 hybridized
                        reactant_num_carbon_neighbors += 1 #we consider it to affect hindrance
                    elif hot_reactant.mol_graph.get_coordination_of_site(neighbor_index) == 1:
                        reactant_num_hydrogens += 1
                        if reactant_num_hydrogens == 3:
                            reaction_methyl_test.append(atom)

        product_num_carbon_neighbors = 0
        for atom in hot_product_atoms: #repeat for products
            product_num_hydrogens = 0
            if hot_product.mol_graph.get_coordination_of_site(atom) == 4:
                neighbor_list = hot_product.mol_graph.get_connected_sites(atom)
                for neighbor in neighbor_list:
                    neighbor_index = neighbor[2]
                    if hot_product.mol_graph.get_coordination_of_site(neighbor_index) == 4:
                        product_num_carbon_neighbors += 1
                    elif hot_product.mol_graph.get_coordination_of_site(neighbor_index) == 1:
                        product_num_hydrogens += 1
                        if product_num_hydrogens == 3:
                            reaction_methyl_test.append(atom)

        if reactant_num_carbon_neighbors >= 3 and product_num_carbon_neighbors >= 3 and len(reaction_methyl_test) < 2: #6 was chosen as the cutoff to prevent tertiary/quaternary carbons from reacting
            return True

        return False


default_reaction_decision_tree = [
    (metal_metal_reaction(), Terminal.DISCARD),
    # redox branch
    (
        is_redox_reaction(),
        [
            (too_many_reactants_or_products(), Terminal.DISCARD),
            (dcharge_too_large(), Terminal.DISCARD),
            (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
            (dG_above_threshold(0.0, "free_energy", 0.0), Terminal.DISCARD),
            (reaction_default_true(), Terminal.KEEP),
        ],
    ),
    (dG_above_threshold(0.0, "solvation_free_energy", 0.0), Terminal.DISCARD),
    # (single_reactant_with_ring_break_two(), Terminal.KEEP),
    # (single_product_with_ring_form_two(), Terminal.KEEP),
    (star_count_diff_above_threshold(6), Terminal.DISCARD),
    (reaction_is_covalent_decomposable(), Terminal.DISCARD),
    (concerted_metal_coordination(), Terminal.DISCARD),
    (concerted_metal_coordination_one_product(), Terminal.DISCARD),
    (concerted_metal_coordination_one_reactant(), Terminal.DISCARD),
    (metal_coordination_passthrough(), Terminal.KEEP),
    (
        fragment_matching_found(),
        [
            (single_reactant_single_product_not_atom_transfer(), Terminal.DISCARD),
            (single_reactant_double_product_ring_close(), Terminal.DISCARD),
            (reaction_default_true(), Terminal.KEEP),
        ],
    ),
    (reaction_default_true(), Terminal.DISCARD),
]

co2_reaction_decision_tree = [
    (
        is_redox_reaction(),
        [
            (too_many_reactants_or_products(), Terminal.DISCARD),
            (dcharge_too_large(), Terminal.DISCARD),
            (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
            (dG_above_threshold(0.0, "free_energy", 0.0), Terminal.DISCARD),
            (reaction_default_true(), Terminal.KEEP),
        ],
    ),
    (dG_above_threshold(0.0, "free_energy", 0.0), Terminal.DISCARD),
    (star_count_diff_above_threshold(6), Terminal.DISCARD),
    (reaction_is_covalent_decomposable(), Terminal.DISCARD),
    (
        fragment_matching_found(),
        [
            (single_reactant_single_product_not_atom_transfer(), Terminal.DISCARD),
            (single_reactant_double_product_ring_close(), Terminal.DISCARD),
            (reaction_default_true(), Terminal.KEEP),
        ],
    ),
    (reaction_default_true(), Terminal.DISCARD),
]


euvl_phase1_reaction_decision_tree_orig = [
    (
        is_redox_reaction(),
        [
            (too_many_reactants_or_products(), Terminal.DISCARD),
            (dcharge_too_large(), Terminal.DISCARD),
            (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
            (add_electron_species(), Terminal.DISCARD),
            (dG_above_threshold(-float("inf"), "free_energy", 0.0), Terminal.KEEP),
            (reaction_default_true(), Terminal.DISCARD),
        ],
    ),
    (reaction_default_true(), Terminal.DISCARD),
]


euvl_phase1_reaction_decision_tree = [
    (
        is_redox_reaction(),
        [
            (too_many_reactants_or_products(), Terminal.DISCARD),
            (dcharge_too_large(), Terminal.DISCARD),
            (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
            (add_electron_species(), Terminal.DISCARD),
            (dG_above_threshold(-float("inf"), "free_energy", 0.0), Terminal.KEEP),
            (reaction_default_true(), Terminal.DISCARD),
        ],
    ),
    (dG_below_threshold(0.0, "free_energy", 0.0), Terminal.DISCARD),
    (
        more_than_one_reactant(), 
        [
            (only_one_product(), Terminal.DISCARD),
            (reactants_are_both_anions_or_both_cations(), Terminal.DISCARD),
            (two_closed_shell_reactants_and_two_open_shell_products(), Terminal.DISCARD),
            (reaction_is_charge_separation(), Terminal.DISCARD),
            (reaction_is_covalent_decomposable(), Terminal.DISCARD),
            (star_count_diff_above_threshold(6), Terminal.DISCARD),
            (compositions_preclude_h_transfer(), Terminal.DISCARD),
            (
                fragment_matching_found(),
                [
                    (not_h_transfer(), Terminal.DISCARD),
                    (h_abstraction_from_closed_shell_reactant(), Terminal.DISCARD),
                    (h_minus_abstraction(), Terminal.DISCARD),
                    (dG_above_threshold(0.0, "free_energy", 0.0, 0.1), Terminal.KEEP),
                    (reaction_default_true(), Terminal.DISCARD),
                ],
            ),
            (reaction_default_true(), Terminal.DISCARD),
        ],
    ),
    (single_reactant_single_product(), Terminal.DISCARD),
    (star_count_diff_above_threshold(4), Terminal.DISCARD),
    (reaction_is_radical_separation(), Terminal.DISCARD),
    (reaction_is_charge_separation(), Terminal.DISCARD),
    (
        fragment_matching_found(),
        [
            (single_reactant_double_product_ring_close(), Terminal.DISCARD),
            (dG_above_threshold(0.0, "free_energy", 0.0), Terminal.KEEP),
            (reaction_default_true(), Terminal.DISCARD),
        ],
    ),
    
    (reaction_default_true(), Terminal.DISCARD),
]

euvl_phase1_reaction_decision_tree_orig = [
    (
        is_redox_reaction(),
        [
            (too_many_reactants_or_products(), Terminal.DISCARD),
            (dcharge_too_large(), Terminal.DISCARD),
            (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
            (add_electron_species(), Terminal.DISCARD),
            (dG_above_threshold(-float("inf"), "free_energy", 0.0), Terminal.KEEP),
            (reaction_default_true(), Terminal.DISCARD),
        ],
    ),
    (reaction_default_true(), Terminal.DISCARD),
]

euvl_phase1_reaction_logging_tree = [
    (
        is_redox_reaction(),
        [
            (too_many_reactants_or_products(), Terminal.DISCARD),
            (dcharge_too_large(), Terminal.DISCARD),
            (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
            (add_electron_species(), Terminal.DISCARD),
            (dG_above_threshold(-float("inf"), "free_energy", 0.0), Terminal.DISCARD),
            (reaction_default_true(), Terminal.DISCARD),
        ],
    ),
    (dG_below_threshold(0.0, "free_energy", 0.0), Terminal.DISCARD),
    (
        more_than_one_reactant(), 
        [
            (only_one_product(), Terminal.DISCARD),
            (reactants_are_both_anions_or_both_cations(), Terminal.DISCARD),
            (two_closed_shell_reactants_and_two_open_shell_products(), Terminal.KEEP),
            (reaction_is_charge_separation(), Terminal.DISCARD),
            (reaction_is_covalent_decomposable(), Terminal.DISCARD),
            (star_count_diff_above_threshold(6), Terminal.DISCARD),
            (compositions_preclude_h_transfer(), Terminal.DISCARD),
            (
                fragment_matching_found(),
                [
                    (not_h_transfer(), Terminal.DISCARD),
                    (h_abstraction_from_closed_shell_reactant(), Terminal.DISCARD),
                    (h_minus_abstraction(), Terminal.DISCARD),
                    (dG_above_threshold(0.0, "free_energy", 0.0, 0.1), Terminal.DISCARD),
                    (reaction_default_true(), Terminal.DISCARD),
                ],
            ),
            (reaction_default_true(), Terminal.DISCARD),
        ],
    ),
    (single_reactant_single_product(), Terminal.DISCARD),
    (star_count_diff_above_threshold(4), Terminal.DISCARD),
    (reaction_is_radical_separation(), Terminal.DISCARD),
    (reaction_is_charge_separation(), Terminal.DISCARD),
    (
        fragment_matching_found(),
        [
            (single_reactant_double_product_ring_close(), Terminal.DISCARD),
            (dG_above_threshold(0.0, "free_energy", 0.0), Terminal.DISCARD),
            (reaction_default_true(), Terminal.DISCARD),
        ],
    ),
    
    (reaction_default_true(), Terminal.DISCARD),
]

euvl_phase2_reaction_decision_tree = [
    (is_redox_reaction(), Terminal.DISCARD),
    (dG_above_threshold(0.0, "free_energy", 0.0), Terminal.DISCARD),
    (reactants_are_both_anions_or_both_cations(), Terminal.DISCARD),
    (reaction_is_charge_transfer(), Terminal.KEEP),
    (reaction_is_covalent_decomposable(), Terminal.DISCARD),
    (star_count_diff_above_threshold(6), Terminal.DISCARD),
    (
        fragment_matching_found(),
        [
            (single_reactant_single_product_not_atom_transfer(), Terminal.DISCARD),
            (single_reactant_double_product_ring_close(), Terminal.DISCARD),
            (reaction_is_hindered(), Terminal.DISCARD),
            (reaction_default_true(), Terminal.KEEP),
        ],
    ),
    (reaction_default_true(), Terminal.DISCARD),
]

euvl_phase2_steric_filter_logging_tree = [
    (is_redox_reaction(), Terminal.DISCARD),
    (dG_above_threshold(0.0, "free_energy", 0.0), Terminal.DISCARD),
    (reactants_are_both_anions_or_both_cations(), Terminal.DISCARD),
    (reaction_is_charge_transfer(), Terminal.DISCARD),
    (reaction_is_covalent_decomposable(), Terminal.DISCARD),
    (star_count_diff_above_threshold(6), Terminal.DISCARD),
    (
        fragment_matching_found(),
        [
            (single_reactant_single_product_not_atom_transfer(), Terminal.DISCARD),
            (single_reactant_double_product_ring_close(), Terminal.DISCARD),
            (reaction_is_hindered(), Terminal.KEEP),
            (reaction_default_true(), Terminal.DISCARD),
        ],
    ),
    (reaction_default_true(), Terminal.DISCARD),
]