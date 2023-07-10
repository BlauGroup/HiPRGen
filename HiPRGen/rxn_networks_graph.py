import json
import os
from pathlib import Path
import copy
from collections import defaultdict
from monty.serialization import dumpfn
from bondnet.data.utils import create_rxn_graph

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

        
        atom_map = rxn['atom_map']
        transformed_atom_map = []
        num_reactants = rxn['number_of_reactants']
        num_products = rxn['number_of_products'] 
        reactants_entry_ids = []
        products_entry_ids = []

        # step 1: Transform atom mapping
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
        assert num_reactants == len(reactants)
        assert num_products == len(products)
        num_tot_atoms = sum([len(i) for i in reactants])
        total_atoms = [i for i in range(num_tot_atoms)]
        

        # step 2: Get total_bonds which is a union of bonds in reactants and products
        reactants_total_bonds = set()
        for k, ind in enumerate(rxn['reactants']):
            mol_reactant = self.mol_entries[ind]
            networkx_graph = mol_reactant.graph
            if len(reactants) <= k: # list index out of range
                break
            reactants_entry_ids.append(self.mol_entries[ind].entry_id)
            for i, j, weight in networkx_graph.edges:
                # reactants_total_bonds.append((k, i, j))
                # print(f"reactants k: {k}")
                # print(f"reactant i: {i}")
                # print(f"reactant j: {j}")
                # print(f"reactants: {reactants}")
                reactants_total_bonds.add(tuple(sorted([reactants[k][i], reactants[k][j]])))
            

        #print(f"reactants_total_bonds: {reactants_total_bonds}")

        products_total_bonds = set()
        for k, ind in enumerate(rxn['products']):
            mol_reactant = self.mol_entries[ind]
            networkx_graph = mol_reactant.graph
            if len(products) <= k: # list index out of range
                break
            products_entry_ids.append(self.mol_entries[ind].entry_id)
            for i, j, weight in networkx_graph.edges:
                # products_total_bonds.append((k, i, j))
                # print(f"products k: {k}")
                # print(f"products i: {i}")
                # print(f"products j: {j}")
                # print(f"products: {products}")
                products_total_bonds.add(tuple(sorted([products[k][i], products[k][j]])))
        
        set_total_bonds = reactants_total_bonds.union(products_total_bonds)
        total_bonds = [[i,j] for i, j in set_total_bonds]
        total_bonds_map = {}
        for ind, bonds in enumerate(total_bonds):
            i, j = bonds
            total_bonds_map[(i,j)] = ind

        
        if rxn['is_redox']:
            assert len(set(reactants_total_bonds)) == len(set(products_total_bonds))
            
        # step 3: Get bond_mapping
        bond_mapping = []
        bonds_in_reactants = [{} for _ in range(num_reactants)]
        for k, ind in enumerate(rxn['reactants']):
            mol_reactant = self.mol_entries[ind]
            networkx_graph = mol_reactant.graph
            if len(reactants) <= k: # list index out of range
                break
            
            for bond_ind, edges in enumerate(networkx_graph.edges):
                i, j, _ = edges
                bonds_in_reactants[k][bond_ind] = total_bonds_map[tuple(sorted([reactants[k][i],reactants[k][j]]))]
        bond_mapping.append(bonds_in_reactants)

        bonds_in_products = [{} for _ in range(num_products)]
        for k, ind in enumerate(rxn['products']):
            mol_reactant = self.mol_entries[ind]
            networkx_graph = mol_reactant.graph
            if len(products) <= k: # list index out of range
                break 
            for bond_ind, edges in enumerate(networkx_graph.edges):
                i, j, _ = edges
                bonds_in_products[k][bond_ind] = total_bonds_map[tuple(sorted([products[k][i],products[k][j]]))]
        bond_mapping.append(bonds_in_products)

        # step 4: get mapping
        mappings = {}
        mappings['bond_map'] = bond_mapping
        mappings['atom_map'] = transformed_atom_map
        mappings['total_bonds'] = total_bonds
        mappings['total_atoms'] = total_atoms
        mappings['num_bonds_total'] = len(total_bonds_map)
        mappings['num_atoms_total'] = len(total_atoms)

        #print(f"mapping: {mappings}")
        print(f"atom_map: {atom_map}")
        print(f"reactants: {reactants}")
        print(f"products: {products}")
        # step 5: Create a reaction graphs and features
        reactants_dgl_graphs  = [self.dgl_mol_dict[entry_i] for entry_i in reactants_entry_ids]
        print(f"reactants_dgl_graphs: {reactants_dgl_graphs}")
        products_dgl_graphs = [self.dgl_mol_dict[entry_i] for entry_i in products_entry_ids]
        print(f"products_dgl_graphs: {products_dgl_graphs}")

        # create has_bonds
        # "has_bonds" is required input to create a reaction graph from BonDNet
        # it's a dictionary e.g) {'reactants': [True, True], 'products': [True]}
        has_bonds = defaultdict(list)
        for _ in range(len(reactants)):
            has_bonds['reactants'].append(True)
        for _ in range(len(products)):
            has_bonds['products'].append(True)

        print(f"has_bonds: {has_bonds}")
        print(f"mappings: {mappings}")
        rxn_graph, features = create_rxn_graph(
                                                reactants = reactants_dgl_graphs,
                                                products = products_dgl_graphs,
                                                mappings = mappings,
                                                has_bonds = has_bonds,
                                                device = None,
                                                ntypes=("global", "atom", "bond"),
                                                ft_name="feat",
                                                reverse=False,
                                            )

        print(f"rxn_graph: {rxn_graph}")
        # print(f"features: {features}")





        # print(f"bond_mapping: {bond_mapping}")
        

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