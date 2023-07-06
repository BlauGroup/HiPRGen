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

        # Transfrom reaction data to the following:
        """mappings = {
            "bond_map": rxn.bond_mapping,
            "atom_map": rxn.atom_mapping,
            "total_bonds": list of list, 
            "total_atoms": list of integer, 
            "num_bonds_total": integer == len(total_bonds),
            "num_atoms_total": integer == len(total_atoms),
        }"""

        # step 1: Transform atom mapping
        atom_map = rxn['atom_map']
        transformed_atom_map = []
        num_reactants = rxn['number_of_reactants']
        num_products = rxn['number_of_products']

        # find the number of atoms for reactant 0
        num_reactant0_atoms = self.mol_entries[rxn['reactants'][0]].num_atoms

        reactants = [{} for _ in range(num_reactants)]
        for ind, atom_i in atom_map.keys():
            reactants[ind][atom_i] = atom_i + ind*num_reactant0_atoms
        transformed_atom_map.append(reactants)

        products = [{} for _ in range(num_products)]
        for r_tuple, p_tuple in atom_map.items():
            prod_ind, p_atom_i = p_tuple
            if r_tuple == p_tuple: 
                products[prod_ind][p_atom_i] = reactants[prod_ind][p_atom_i]
            else:
                react_ind, r_atom_i = r_tuple
                products[prod_ind][p_atom_i] = reactants[react_ind][r_atom_i]
        transformed_atom_map.append(products)
        # print(f"products: {products}")
        # print(f"reactants: {reactants}")
        # print(f"transformed_atom_map: {transformed_atom_map}")
        # check the conservation of mass in a reaction
        assert sum([len(i) for i in reactants]) == sum([len(i) for i in products])

        # step 2: Get total_bonds
        reactants_total_bonds = []
        for k, ind in enumerate(rxn['reactants']):
            mol_reactant = self.mol_entries[ind]
            networkx_graph = mol_reactant.graph
            for i, j, weight in networkx_graph.edges:
                reactants_total_bonds.append((k, i, j))
                #reactants_total_bonds.append((reactants[k][i], reactants[k][j]))

        #print(f"reactants_total_bonds: {reactants_total_bonds}")
        len_reactants_total_bonds = len(reactants_total_bonds)

        products_total_bonds = []
        for k, ind in enumerate(rxn['products']):
            mol_reactant = self.mol_entries[ind]
            networkx_graph = mol_reactant.graph
            for i, j, weight in networkx_graph.edges:
                products_total_bonds.append((k, i, j))
                #products_total_bonds.append((products[k][i], products[k][j]))

        # 
        if rxn['is_redox']:
            assert len(set(reactants_total_bonds)) == len(set(products_total_bonds))
        else:
            print(f"reactants_total_bonds: {reactants_total_bonds}")
            print(f"products_total_bonds: {products_total_bonds}")
            print(f'transformed_atom_map: {transformed_atom_map}')
            print(f"reactant bonds broken: {rxn['reactant_bonds_broken']}")
            print(f"product bonds broken: {rxn['product_bonds_broken']}")

        #print(f"products_total_bonds: {products_total_bonds}")

        #assert len(products_total_bonds) == len_reactants_total_bonds
        # reactant_bonds_broken = rxn["reactant_bonds_broken"]
        # product_bonds_broken = rxn["product_bonds_broken"]
        # print(f"reaction bonds broken: {reactant_bonds_broken}")
        # print(f"reaction bonds broken: {product_bonds_broken}")


        #bonds_intersection = reactants_total_bonds.intersection(products_total_bonds)

        # assert len(bonds_intersection) == len_reactants_total_bonds-1
        # assert len(bonds_intersection) == len(products_total_bonds)-1

        #total_bonds = reactants_total_bonds.union(products_total_bonds)
    
        # print(f"total bonds: {total_bonds}")

            


        # step 3: Get bond_map
        

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