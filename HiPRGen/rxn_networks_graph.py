import torch
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
        
        # initialize JSON data
        self.data = {}

    def create_rxn_networks_graph(self, rxn, rxn_id):

        # Transfrom reaction data to the following:
        """mappings = {
            "bond_map": a list of lists, 
                        e.g.) For a reaction that has two reactants and two products, it has [[{}, {}], [{}, {}]] format.
                        The first inner list includes dictionaries of bond map for reactants.
                        The second inner list includes dictionaries of bond map for products.
                        A key of the dictionary represents a local bond index and the corresponding value represents
                        a global bond index from "total_bonds".
            "atom_map": a list of lists,
                        e.g.) For a reaction that has two reactants and two products, it has [[{}, {}], [{}, {}]] format.
                        The first inner list includes dictionaries of atom map for reactants.
                        The second inner list includes dictionaries of atom map for products.
                        A key of the dictionary represetns a local atom index and the corresponding value represents 
                        a global atom index.
            "total_bonds": list of lists whose length is 2. It's an union of bonds in reactants and products
                        e.g.) [[0,1], [0,2], [0,3], [0, 4]]
            "total_atoms": list of integer, 
            "num_bonds_total": integer == len(total_bonds),
            "num_atoms_total": integer == len(total_atoms),
        }"""

        
        atom_map = rxn['atom_map']
        num_reactants = rxn['number_of_reactants']
        num_products = rxn['number_of_products'] 
        transformed_atom_map = []
        
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

        # Get "total_atoms"
        num_tot_atoms = sum([len(i) for i in reactants])
        total_atoms = [i for i in range(num_tot_atoms)]
        

        # step 2: Get total_bonds which is a union of bonds in reactants and products
        def find_total_bonds(rxn, species, reactants, products):
            species_entry_ids = []
            species_total_bonds = set()
            if species == 'reactants':
                temp_species = reactants
            else:
                temp_species = products
            for k, ind in enumerate(rxn[species]):
                mol_reactant = self.mol_entries[ind]
                networkx_graph = mol_reactant.graph
                # This is needed because there is a case where num_reactants != len(rxn['reactants']) or num_products != len(rxn['products'])
                if len(temp_species) <= k: 
                    break
                species_entry_ids.append(self.mol_entries[ind].entry_id)
                for i, j, weight in networkx_graph.edges:
                    species_total_bonds.add(tuple(sorted([temp_species[k][i], temp_species[k][j]])))
            return species_total_bonds, species_entry_ids
        
        reactants_total_bonds, reactants_entry_ids = find_total_bonds(rxn, 'reactants', reactants, products)
        products_total_bonds, products_entry_ids = find_total_bonds(rxn, 'products', reactants, products)

        # reactants_total_bonds = set()
        # for k, ind in enumerate(rxn['reactants']):
        #     mol_reactant = self.mol_entries[ind]
        #     networkx_graph = mol_reactant.graph
        #     # This is needed because there is a case where num_reactants != len(rxn['reactants'])
        #     if len(reactants) <= k: 
        #         break
        #     reactants_entry_ids.append(self.mol_entries[ind].entry_id)
        #     for i, j, weight in networkx_graph.edges:
        #         reactants_total_bonds.add(tuple(sorted([reactants[k][i], reactants[k][j]])))
            

        # products_total_bonds = set()
        # for k, ind in enumerate(rxn['products']):
        #     mol_reactant = self.mol_entries[ind]
        #     networkx_graph = mol_reactant.graph
        #     if len(products) <= k: # list index out of range
        #         break
        #     products_entry_ids.append(self.mol_entries[ind].entry_id)
        #     for i, j, weight in networkx_graph.edges:
        #         products_total_bonds.add(tuple(sorted([products[k][i], products[k][j]])))
        
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
        # print(f"atom_map: {atom_map}")
        # print(f"reactants: {reactants}")
        # print(f"products: {products}")
        # step 5: Create a reaction graphs and features
        reactants_dgl_graphs  = [self.dgl_mol_dict[entry_i] for entry_i in reactants_entry_ids]
        # print(f"reactants_dgl_graphs: {reactants_dgl_graphs}")
        products_dgl_graphs = [self.dgl_mol_dict[entry_i] for entry_i in products_entry_ids]
        # print(f"products_dgl_graphs: {products_dgl_graphs}")

        # create has_bonds
        # "has_bonds" is required input to create a reaction graph from BonDNet
        # it's a dictionary e.g) {'reactants': [True, True], 'products': [True]}
        has_bonds = defaultdict(list)
        for _ in range(len(reactants)):
            has_bonds['reactants'].append(True)
        for _ in range(len(products)):
            has_bonds['products'].append(True)

        # print(f"has_bonds: {has_bonds}")
        # print(f"mappings: {mappings}")
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

        #print(f"rxn_graph: {rxn_graph}")
        # print(f"features: {features}")

        # step 5: update reaction features to the reaction graph
        for nt, ft in features.items():
            rxn_graph.nodes[nt].data.update({'ft': ft})

        self.data[rxn_id] = {}
        self.data[rxn_id]['rxn_graph'] = str(rxn_graph)
        self.data[rxn_id]['value'] = str(torch.tensor([rxn['dG']]))

        

    def save_data(self):
        #self.f.write(json.dumps(self.data))
        dumpfn(self.data, self.report_file_path)