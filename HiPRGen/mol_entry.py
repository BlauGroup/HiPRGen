import copy
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
from monty.json import MSONable
from pymatgen.analysis.graphs import MoleculeGraph, MolGraphSplitError
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender
from pymatgen.core.structure import Molecule
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from HiPRGen.constants import ROOM_TEMP
from itertools import permutations, product


metals = frozenset(["Li", "Na", "K", "Mg", "Ca", "Zn", "Al"])
m_formulas = frozenset([m + "1" for m in metals])

solvation_correction = {
    "Li" : -0.68,
    }

coordination_radius = {
    "Li" : 2.4,
    }

max_number_of_coordination_bonds = {
    "Li" : 4,
    }


class Fragment:
    def __init__(
            self,
            fragment_hash,
            atom_ids,
            neighborhood_hashes,
            graph,
            hot_atoms
    ):

        self.fragment_hash = fragment_hash
        self.atom_ids = atom_ids
        self.neighborhood_hashes = neighborhood_hashes
        self.graph = graph
        self.hot_atoms = hot_atoms


def sym_iterator(n):
    return permutations(range(n), r=n)


def find_fragment_atom_mappings(fragment_1, fragment_2):
    groups_by_hash = {}

    for left_index in fragment_1.atom_ids:

        neighborhood_hash = fragment_1.neighborhood_hashes[left_index]
        if neighborhood_hash not in groups_by_hash:
            groups_by_hash[neighborhood_hash] = ([],[])

        groups_by_hash[neighborhood_hash][0].append(left_index)


    for right_index in fragment_2.atom_ids:

        neighborhood_hash = fragment_2.neighborhood_hashes[right_index]
        if neighborhood_hash not in groups_by_hash:
            groups_by_hash[neighborhood_hash] = ([],[])

        groups_by_hash[neighborhood_hash][1].append(right_index)

    groups = list(groups_by_hash.values())

    product_sym_iterator = product(*[
        sym_iterator(len(p[0])) for
        p in groups ])

    mappings = []

    for product_perm in product_sym_iterator:
        mapping = {}
        for perm, vals in zip(product_perm, groups):
            for i, j in enumerate(perm):
                mapping[vals[0][i]] = vals[1][j]

        isomorphism = True
        for edge in fragment_1.graph.edges:
            u = mapping[edge[0]]
            v = mapping[edge[1]]
            if not fragment_2.graph.has_edge(u,v):
                isomorphism = False
                break


        if isomorphism:
            mappings.append(mapping)

    return mappings

def find_hot_atom_preserving_fragment_map(
        fragment_1,
        fragment_2,
        mappings):

    for mapping in mappings:
        hot_atom_preserving = False
        for hot_atom in fragment_1.hot_atoms:
            if mapping[hot_atom] in fragment_2.hot_atoms:
                hot_atom_preserving = True
                break

        if hot_atom_preserving:
            return mapping

    return None



class FragmentComplex:

    def __init__(
            self,
            number_of_fragments,
            number_of_bonds_broken,
            fragments):

        self.number_of_fragments = number_of_fragments
        self.number_of_bonds_broken = number_of_bonds_broken
        self.fragments = fragments



class MoleculeEntry(MSONable):
    """
    A molecule entry class to provide easy access to Molecule properties.

    Args:
        molecule: Molecule of interest.
        energy: Electronic energy of the molecule in Hartree.
        enthalpy: Enthalpy of the molecule (kcal/mol). Defaults to None.
        entropy: Entropy of the molecule (cal/mol.K). Defaults to None.
        entry_id: An optional id to uniquely identify the entry.
        mol_graph: MoleculeGraph of the molecule.
    """

    def __init__(
        self,
        molecule: Molecule,
        energy: float,
        enthalpy: Optional[float] = None,
        entropy: Optional[float] = None,
        entry_id: Optional[Any] = None,
        mol_graph: Optional[MoleculeGraph] = None,
        partial_charges_resp: Optional[list] = None,
        partial_charges_mulliken: Optional[list] = None
    ):
        self.energy = energy
        self.enthalpy = enthalpy
        self.entropy = entropy


        self.ind = None
        self.entry_id = entry_id

        self.star_hashes = {}
        self.fragment_data = []


        if not mol_graph:
            mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())
            self.mol_graph = metal_edge_extender(mol_graph)
        else:
            self.mol_graph = mol_graph

        self.partial_charges_resp = partial_charges_resp
        self.partial_charges_mulliken = partial_charges_mulliken
        self.molecule = self.mol_graph.molecule
        self.graph = self.mol_graph.graph.to_undirected()
        self.species = [str(s) for s in self.molecule.species]

        self.m_inds = [
            i for i, x in enumerate(self.species) if x in metals
        ]


        self.covalent_graph = copy.deepcopy(self.graph)
        self.covalent_graph.remove_nodes_from(self.m_inds)


        self.total_hash = weisfeiler_lehman_graph_hash(
                self.graph,
                node_attr='specie')

        self.covalent_hash = weisfeiler_lehman_graph_hash(
                self.covalent_graph,
                node_attr='specie')


        self.formula = self.molecule.composition.alphabetical_formula
        self.charge = self.molecule.charge
        self.num_atoms = len(self.molecule)

        self.atom_locations = [
            site.coords for site in self.molecule]


        self.number_of_coordination_bonds = 0
        self.free_energy = self.get_free_energy()
        self.solvation_free_energy = self.get_solvation_free_energy()

        self.non_metal_atoms = [
            i for i in range(self.num_atoms)
            if self.species[i] not in metals]




    @classmethod
    def from_dataset_entry(
        cls,
        doc: Dict,
        use_thermo: str = "raw",
    ):
        """
        Initialize a MoleculeEntry from a document in the LIBE (Lithium-Ion
        Battery Electrolyte) or MADEIRA (MAgnesium Dataset of Electrolyte and
        Interphase ReAgents) datasets.

        Args:
            doc: Dictionary representing an entry from LIBE or MADEIRA
            use_thermo: One of "raw" (meaning raw, uncorrected thermo data will
                be used), "rrho_shifted" (meaning that a slightly modified
                Rigid-Rotor Harmonic Oscillator approximation will be used -
                see Ribiero et al., J. Phys. Chem. B 2011, 115, 14556-14562), or
                "qrrho" (meaning that Grimme's Quasi-Rigid Rotor Harmonic
                Oscillator - see Grimme, Chem. Eur. J. 2012, 18, 9955-9964) will
                be used.
        """

        thermo = use_thermo.lower()

        if thermo not in ["raw", "rrho_shifted", "qrrho"]:
            raise ValueError(
                "Only allowed values for use_thermo are 'raw', 'rrho_shifted', "
                "and 'qrrho'!"
            )
        try:
            if isinstance(doc["molecule"], Molecule):
                molecule = doc["molecule"]
            else:
                molecule = Molecule.from_dict(doc["molecule"])  # type: ignore

            if (
                thermo == "rrho_shifted"
                and doc["thermo"]["shifted_rrho_eV"] is not None
            ):
                energy = (
                    doc["thermo"]["shifted_rrho_eV"]["electronic_energy"] * 0.0367493
                )
                enthalpy = doc["thermo"]["shifted_rrho_eV"]["total_enthalpy"] * 23.061
                entropy = doc["thermo"]["shifted_rrho_eV"]["total_entropy"] * 23061
            elif thermo == "qrrho" and doc["thermo"]["quasi_rrho_eV"] is not None:
                energy = doc["thermo"]["quasi_rrho_eV"]["electronic_energy"] * 0.0367493
                enthalpy = doc["thermo"]["quasi_rrho_eV"]["total_enthalpy"] * 23.061
                entropy = doc["thermo"]["quasi_rrho_eV"]["total_entropy"] * 23061
            else:
                energy = doc["thermo"]["raw"]["electronic_energy_Ha"]
                enthalpy = doc["thermo"]["raw"]["total_enthalpy_kcal/mol"]
                entropy = doc["thermo"]["raw"]["total_entropy_cal/molK"]

            entry_id = doc["molecule_id"]

            if isinstance(doc["molecule_graph"], MoleculeGraph):
                mol_graph = doc["molecule_graph"]
            else:
                mol_graph = MoleculeGraph.from_dict(doc["molecule_graph"])

            partial_charges_resp = doc['partial_charges']['resp']
            partial_charges_mulliken = doc['partial_charges']['mulliken']
        except KeyError as e:
            raise MoleculeEntryError(
                "Unable to construct molecule entry from molecule document; missing "
                f"attribute {e} in `doc`."
            )



        return cls(
            molecule=molecule,
            energy=energy,
            enthalpy=enthalpy,
            entropy=entropy,
            entry_id=entry_id,
            mol_graph=mol_graph,
            partial_charges_resp=partial_charges_resp,
            partial_charges_mulliken=partial_charges_mulliken
        )



    def get_free_energy(self, temperature: float = ROOM_TEMP) -> Optional[float]:
        """
        Get the free energy at the give temperature.
        """
        if self.enthalpy is not None and self.entropy is not None:
            # TODO: fix these hard coded vals
            return (
                self.energy * 27.21139
                + 0.0433641 * self.enthalpy
                - temperature * self.entropy * 0.0000433641
            )
        else:
            return None


    def get_solvation_free_energy(self):
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

        # this method should only be called once, but just to be safe,
        # reset the coordination bond count
        self.number_of_coordination_bonds = 0

        correction = 0.0

        for i in self.m_inds:

            species = self.species[i]
            coordination_partners = []
            radius = coordination_radius[species]

            for j in range(self.num_atoms):
                if j != i:
                    displacement_vector = (
                        self.atom_locations[j] -
                        self.atom_locations[i])

                    if (np.inner(displacement_vector, displacement_vector)
                        < radius ** 2 and (
                            self.partial_charges_resp[j] < 0 or
                            self.partial_charges_mulliken[j] < 0)):
                        coordination_partners.append(j)


            number_of_coordination_bonds = len(coordination_partners)
            self.number_of_coordination_bonds += number_of_coordination_bonds
            correction += solvation_correction[species] * (
                max_number_of_coordination_bonds[species] -
                number_of_coordination_bonds)

        return correction + self.free_energy


    def __repr__(self):

        output = [
            f"MoleculeEntry {self.entry_id} - {self.formula}",
            f"Total charge = {self.charge}",
        ]

        energies = [
            ("Energy", "Hartree", self.energy),
            ("Enthalpy", "kcal/mol", self.enthalpy),
            ("Entropy", "cal/mol.K", self.entropy),
            ("Free Energy (298.15 K)", "eV", self.get_free_energy()),
        ]
        for name, unit, value in energies:
            if value is None:
                output.append(f"{name} = {value} {unit}")
            else:
                output.append(f"{name} = {value:.4f} {unit}")

        if self.ind:
            output.append("index: {}".format(self.ind))

        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if type(self) == type(other):
            return str(self) == str(other)
        else:
            return False




class MoleculeEntryError(Exception):
    def __init__(self, message):
        super().__init__(message)
        self.message = message
