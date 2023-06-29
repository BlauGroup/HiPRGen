import json
import os
from pathlib import Path
import copy
from monty.serialization import dumpfn

class rxn_networks_graph:
    def __init__(
        self,
        mol_entries,
        dgl_molecules_dict,
        report_file_path
    ):
        self.mol_entries = mol_entries
        self.dgl_mol_dict = dgl_molecules_dict
        self.report_file_path = report_file_path
        
        self.keys = ['atom_map', 'bond_map', 'total_bonds', 'total_atoms', 'num_bonds_total', 'num_atoms_total']
        
        # initialize JSON data
        self.data = {}
        for k in self.keys:
            self.data[k] = {}

    def create_rxn_networks_graph(self, rxn, rxn_id):

        # step 1: transfrom reaction data to the following:
        """mappings = {
            "bond_map": rxn.bond_mapping,
            "atom_map": rxn.atom_mapping,
            "total_bonds": rxn.total_bonds, 
            "total_atoms": rxn.total_atoms, 
            "num_bonds_total": rxn.num_bonds_total,
            "num_atoms_total": rxn.num_atoms_total,
        }"""
        self.data['atom_map'][rxn_id] = str(rxn['atom_map'])
        print("printing dgl_mol_dict")
        print(self.dgl_mol_dict)



    # def insert_data(self, rxn, rxn_id):
        # should be updated
        # assume only one reactant for now 
        # mol_reactant = self.mol_entries[rxn["reactants"][0]] # MolecularEntry object
        # mol_product = self.mol_entries[rxn["products"][0]]
        # print(mol_reactant.entry_id)
        # graph_hash = mol_reactant.entry_id.split('-')[0]
        # self.data['_id'][rxn_id] = graph_hash
        # self.data['builder_meta'][rxn_id] = None
        # self.data['charge'][rxn_id] = mol_reactant.charge
        # self.data['spin_multiplicity'][rxn_id] = mol_reactant.spin_multiplicity
        # self.data['natoms'][rxn_id] = mol_reactant.num_atoms
        # self.data['elements'][rxn_id] = list(set(mol_reactant.species))
        # self.data['nelements'][rxn_id] = len(self.data['elements'][rxn_id])
        # self.data['nelectrons'][rxn_id] = mol_reactant.molecule.nelectrons
        # self.data['composition'][rxn_id] = mol_reactant.molecule.composition.as_dict()
        # self.data['composition_reduced'][rxn_id] = mol_reactant.molecule.composition.to_reduced_dict
        # self.data['formula_alphabetical'][rxn_id] = mol_reactant.formula
        # self.data['formula_pretty'][rxn_id] = mol_reactant.molecule.composition.to_pretty_string()
        # self.data['formula_anonymous'][rxn_id] = mol_reactant.molecule.composition.anonymized_formula
        # self.data['chemsys'][rxn_id] = mol_reactant.molecule.composition.chemical_system
        # self.data['symmetry'][rxn_id] = None
        # self.data['reactant_structure'][rxn_id] = mol_reactant.molecule.as_dict()
        # self.data['reactant_molecule_graph'][rxn_id] = mol_reactant.mol_graph
    

    def save_data(self):
        #self.f.write(json.dumps(self.data))
        dumpfn(self.data, self.report_file_path)