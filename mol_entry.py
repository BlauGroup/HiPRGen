# coding: utf-8
# Copyright (c) MR.Net Development Team.
# Distributed under the terms of the MIT License.

import copy
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
from monty.json import MSONable
from pymatgen.analysis.graphs import MoleculeGraph, MolGraphSplitError
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender
from pymatgen.core.structure import Molecule

from constants import ROOM_TEMP

__author__ = "Sam Blau, Mingjian Wen"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "Aug 1, 2019"


class MoleculeEntry(MSONable):
    """
    A molecule entry class to provide easy access to Molecule properties.

    Args:
        molecule: Molecule of interest.
        energy: Electronic energy of the molecule in Hartree.
        correction: A correction to be applied to the energy.
            This is used to modify the energy for certain analyses.
            Defaults to 0.0.
        enthalpy: Enthalpy of the molecule (kcal/mol). Defaults to None.
        entropy: Entropy of the molecule (cal/mol.K). Defaults to None.
        parameters: An optional dict of parameters associated with
            the molecule. Defaults to None.
        entry_id: An optional id to uniquely identify the entry.
        attribute: Optional attribute of the entry. This can be used to
            specify that the entry is a newly found compound, or to specify
            a particular label for the entry, or else ... Used for further
            analysis and plotting purposes. An attribute can be anything
            but must be MSONable.
        mol_graph: MoleculeGraph of the molecule.
    """

    def __init__(
        self,
        molecule: Molecule,
        energy: float,
        correction: float = 0.0,
        enthalpy: Optional[float] = None,
        entropy: Optional[float] = None,
        parameters: Optional[Dict] = None,
        entry_id: Optional[Any] = None,
        attribute=None,
        mol_graph: Optional[MoleculeGraph] = None,
    ):
        self.uncorrected_energy = energy
        self.correction = correction
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.parameters = parameters if parameters else {}
        self.entry_id = entry_id
        self.attribute = attribute

        if not mol_graph:
            mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())
            self.mol_graph = metal_edge_extender(mol_graph)
        else:
            self.mol_graph = mol_graph

    @classmethod
    def from_molecule_document(
        cls,
        mol_doc: Dict,
        correction: float = 0.0,
        parameters: Optional[Dict] = None,
        attribute=None,
    ):
        """
        Initialize a MoleculeEntry from a molecule document.

        Args:
            mol_doc: MongoDB molecule document (nested dictionary) that contains the
                molecule information.
            correction: A correction to be applied to the energy. This is used to modify
                the energy for certain analyses. Defaults to 0.0.
            parameters: An optional dict of parameters associated with
                the molecule. Defaults to None.
            attribute: Optional attribute of the entry. This can be used to
                specify that the entry is a newly found compound, or to specify
                a particular label for the entry, or else ... Used for further
                analysis and plotting purposes. An attribute can be anything
                but must be MSONable.
        """
        try:
            if isinstance(mol_doc["molecule"], Molecule):
                molecule = mol_doc["molecule"]
            else:
                molecule = Molecule.from_dict(mol_doc["molecule"])  # type: ignore
            energy = mol_doc["energy_Ha"]
            enthalpy = mol_doc["enthalpy_kcal/mol"]
            entropy = mol_doc["entropy_cal/molK"]
            entry_id = mol_doc["task_id"]
        except KeyError as e:
            raise MoleculeEntryError(
                "Unable to construct molecule entry from molecule document; missing "
                f"attribute {e} in `mol_doc`."
            )

        if "mol_graph" in mol_doc:
            if isinstance(mol_doc["mol_graph"], MoleculeGraph):
                mol_graph = mol_doc["mol_graph"]
            else:
                mol_graph = MoleculeGraph.from_dict(mol_doc["mol_graph"])
        else:
            mol_graph = None

        return cls(
            molecule=molecule,
            energy=energy,
            correction=correction,
            enthalpy=enthalpy,
            entropy=entropy,
            parameters=parameters,
            entry_id=entry_id,
            attribute=attribute,
            mol_graph=mol_graph,
        )

    @classmethod
    def from_dataset_entry(
        cls,
        doc: Dict,
        use_thermo: str = "raw",
        parameters: Optional[Dict] = None,
        attribute=None,
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
            parameters: An optional dict of parameters associated with
                the molecule. Defaults to None.
            attribute: Optional attribute of the entry. This can be used to
                specify that the entry is a newly found compound, or to specify
                a particular label for the entry, or else ... Used for further
                analysis and plotting purposes. An attribute can be anything
                but must be MSONable.
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
            parameters=parameters,
            entry_id=entry_id,
            attribute=attribute,
            mol_graph=mol_graph,
        )

    @property
    def molecule(self):
        return self.mol_graph.molecule

    @property
    def graph(self) -> nx.MultiDiGraph:
        return self.mol_graph.graph

    @property
    def energy(self) -> float:
        return self.uncorrected_energy + self.correction

    @property
    def formula(self) -> str:
        return self.molecule.composition.alphabetical_formula

    @property
    def charge(self) -> float:
        return self.molecule.charge

    @property
    def species(self) -> List[str]:
        return [str(s) for s in self.molecule.species]

    @property
    def bonds(self) -> List[Tuple[int, int]]:
        return [(int(sorted(e)[0]), int(sorted(e)[1])) for e in self.graph.edges()]

    @property
    def num_atoms(self) -> int:
        return len(self.molecule)

    @property
    def num_bonds(self) -> int:
        return len(self.bonds)

    @property
    def coords(self) -> np.ndarray:
        return self.molecule.cart_coords

    def get_free_energy(self, temperature: float = ROOM_TEMP) -> Optional[float]:
        """
        Get the free energy at the give temperature.
        """
        if self.enthalpy is not None and self.entropy is not None:
            return (
                self.energy * 27.21139
                + 0.0433641 * self.enthalpy
                - temperature * self.entropy * 0.0000433641
            )
        else:
            return None

    def get_fragments(self) -> Optional[Dict[Tuple[Any, Any], List[MoleculeGraph]]]:
        """
        Get the fragments of the molecule by breaking all its bonds.

        Returns:
            Fragments dict {(atom1, atom2): [fragments]}, where
                the key `(atom1, atom2)` specifies the broken bond indexed by the two
                atoms forming the bond, and the value `[fragments]` is a list of
                fragments obtained by breaking the bond. This list can have either one
                element (ring-opening A->B) or two elements (not ring-opening A->B+C).
                The dictionary is empty if the molecule has no bonds (e.g. Li+).
        """

        fragments = {}
        if self.mol_graph:
            for edge in self.bonds:
                try:
                    frags = self.mol_graph.split_molecule_subgraphs(
                        [edge], allow_reverse=True, alterations=None
                    )
                    fragments[edge] = frags

                except MolGraphSplitError:
                    # cannot split (ring-opening editing)
                    frag = copy.deepcopy(self.mol_graph)
                    idx1, idx2 = edge
                    frag.break_edge(idx1, idx2, allow_reverse=True)
                    fragments[edge] = [frag]

            return fragments
        else:
            return None

    def get_isomorphic_bonds(
        self, fragments: Optional[Dict[Tuple[int, int], List[MoleculeGraph]]] = None
    ) -> Optional[List[List[Tuple[int, int]]]]:
        """
        Find isomorphic bonds in the molecule.

        Isomorphic bonds are defined as bonds that when breaking them separately,
        the same fragments (in terms of graph connectivity) are obtained.

        For example, for molecule:

             b0      b1
        H(1)----C(0)----H(2)
            b2 /   | b3
            O(3)---O(4)
                b4

        (notation: number after b is bond index, number in `()` is atom index)

        bond 0 is isomorphic to bond 1, and bond 2 is isomorphic to bond 3.

        Args:
            fragments: a dictionary of fragments obtained by breaking all bonds in the
                molecule, can be obtained by `self.get_fragments()`. If `None`,
                will generate the fragments automatically.

        Returns:
            Isomorphic bonds specified by atom indexes: [[(atom1, atom2)]].
            Each inner list contains bonds (each bond is specified by the indexes of the
            atoms forming the bond `(atom1, atom2)`) that are isomorphic to each other.
            Note, bond not isomorphic to any other bond is included as a group by itself.
            For example, for the above shown molecule, this function returns:
            [[(0,1), (0,2)], [(0,3), (0,4)], [(3,4)]]
        """

        fragments = self.get_fragments() if fragments is None else fragments

        if fragments:

            iso_bonds = []  # type: List[List[Tuple[int, int]]]

            for current_bond, current_frags in fragments.items():
                for group in iso_bonds:

                    # compare to the first element in a group to determine whether they are
                    # isomorphic to each other
                    existing_bond = group[0]
                    exsiting_frags = fragments[existing_bond]

                    # one fragments (ring-opening like fragments)
                    if len(current_frags) == len(exsiting_frags) == 1:
                        if current_frags[0].isomorphic_to(exsiting_frags[0]):
                            group.append(current_bond)
                            break

                    # two fragments
                    elif len(current_frags) == len(exsiting_frags) == 2:
                        if (
                            current_frags[0].isomorphic_to(exsiting_frags[0])
                            and current_frags[1].isomorphic_to(exsiting_frags[1])
                        ) or (
                            current_frags[0].isomorphic_to(exsiting_frags[1])
                            and current_frags[1].isomorphic_to(exsiting_frags[0])
                        ):
                            group.append(current_bond)
                            break

                # current_bond not in any group, create a new group
                else:
                    iso_bonds.append([current_bond])

            return iso_bonds
        else:
            return None

    def __repr__(self):

        output = [
            f"MoleculeEntry {self.entry_id} - {self.formula}",
            f"Number of bonds = {self.num_bonds}",
            f"Total charge = {self.charge}",
        ]

        energies = [
            ("Energy", "Hartree", self.uncorrected_energy),
            ("Correction", "Hartree", self.correction),
            ("Enthalpy", "kcal/mol", self.enthalpy),
            ("Entropy", "cal/mol.K", self.entropy),
            ("Free Energy (298.15 K)", "eV", self.get_free_energy()),
        ]
        for name, unit, value in energies:
            if value is None:
                output.append(f"{name} = {value} {unit}")
            else:
                output.append(f"{name} = {value:.4f} {unit}")

        if self.parameters:
            output.append("Parameters:")
            for k, v in self.parameters.items():
                output.append("{} = {}".format(k, v))

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
