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


metals = frozenset(["Li", "Na", "K", "Mg", "Ca", "Zn", "Al"])
m_formulas = frozenset([m + "1" for m in metals])

solvation_correction = {
    "Li" : -0.68,
    }

max_number_of_bonds = {
    "Li" : 4,
    }

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
    ):
        self.energy = energy
        self.enthalpy = enthalpy
        self.entropy = entropy


        self.ind = None

        self.star_hashes = {}
        self.fragment_hashes = []



        self.entry_id = entry_id

        if not mol_graph:
            mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())
            self.mol_graph = metal_edge_extender(mol_graph)
        else:
            self.mol_graph = mol_graph

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
        self.bonds = [(int(sorted(e)[0]), int(sorted(e)[1])) for e in self.graph.edges()]
        self.num_atoms = len(self.molecule)
        self.num_bonds = len(self.bonds)

        self.free_energy = self.get_free_energy()
        self.solvation_free_energy = self.get_solvation_free_energy()



    @classmethod
    def from_molecule_document(
        cls,
        mol_doc: Dict,
    ):
        """
        Initialize a MoleculeEntry from a molecule document.

        Args:
            mol_doc: MongoDB molecule document (nested dictionary) that contains the
                molecule information.
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
            enthalpy=enthalpy,
            entropy=entropy,
            entry_id=entry_id,
            mol_graph=mol_graph,
        )

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

        correction = 0.0

        for i in self.m_inds:
            number_of_bonds = len([bond for bond in self.bonds if i in bond])
            species = self.species[i]
            correction += solvation_correction[species] * (
                max_number_of_bonds[species] - number_of_bonds)

        return correction + self.free_energy


    def __repr__(self):

        output = [
            f"MoleculeEntry {self.entry_id} - {self.formula}",
            f"Number of bonds = {self.num_bonds}",
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
