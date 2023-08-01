import torch
import math
import os
from pathlib import Path
import copy
from collections import defaultdict
from monty.serialization import dumpfn
from bondnet.data.utils import create_rxn_graph
from HiPRGen.lmdb_dataset import LmdbDataset
import lmdb
import tqdm
import pickle
from HiPRGen.lmdb_dataset import write_to_lmdb

class rxn_networks_graph:
    def __init__(
        self,
        mol_entries,
        dgl_molecules_dict,
        grapher_features,
        report_file_path
    ):
        self.mol_entries = mol_entries
        self.dgl_mol_dict = dgl_molecules_dict
        self.grapher_features = grapher_features
        self.report_file_path = report_file_path
        
        # initialize data
        self.data = {} 
        

    def create_rxn_networks_graph(self, rxn, rxn_id):

        """
        To create a reaction graph (or a reaction networks graph) using a function from BonDNet,
        we need to specify four inputs:
        - reactants: a list of dgl graphs of reactants,
        - products: a list of dgl graphs of products,
        - mappings: detailed explanation in the below,
        - has_bonds: a dictionary with a form of {'reactants': [True, True], 'products': [True]} as an example of two reactants and one product.
        
        mappings = {
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
            "total_bonds": list of lists. It's an union of bonds in reactants and products
                        e.g.) [[0,1], [0,2], [0,3], [0, 4]]
            "total_atoms": list of integer, 
            "num_bonds_total": integer == len(total_bonds),
            "num_atoms_total": integer == len(total_atoms),
            }

        The goal of this function is to create a reaction graph with a reaction filtered from HiPRGen, specifically, reaction_filter.py
        """

        
        atom_map = rxn['atom_map']
        num_reactants = rxn['number_of_reactants']
        num_products = rxn['number_of_products'] 
        transformed_atom_map = []

        # step 1: Transform atom map to "atom_map" format in mappings

        # find the number of atoms for reactant_0
        num_reactant0_atoms = self.mol_entries[rxn['reactants'][0]].num_atoms

        # transform atom map for reactants
        reactants = [{} for _ in range(num_reactants)]
        for ind, atom_i in atom_map.keys():
            reactants[ind][atom_i] = atom_i + ind*num_reactant0_atoms
        transformed_atom_map.append(reactants)

        # transform atom map for products
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

        # Find "total_atoms" in mapping
        num_tot_atoms = sum([len(i) for i in reactants])
        total_atoms = [i for i in range(num_tot_atoms)]
        

        # step 2: Get total_bonds which are a union of bonds in reactants and products

        # define a function to find total bonds in reactants or products
        def find_total_bonds(rxn, species, reactants, products):
            """ Goal: find total bonds in reactants or products and fetch entry_ids of reactants or products 
                Inputs:
                - rxn: a reaction from reaction_filter.py from HiPRGen
                - species: (str) 'reactants' or 'products'
                - reactants: a list of dictionaries of transformed atom map for reactants
                - products: a list of dictionaries of transformed atom map for products
                Outputs:
                - species_entry_ids: a list of mol.entry_ids of reactants or products
                - species_total_bonds:  a set of tuples. Each tuple represents a bond between two atoms noted with global indexes
            """
            species_entry_ids = []
            species_total_bonds = set()
            original_bonds = []
            if species == 'reactants':
                temp_species = reactants
            else:
                temp_species = products
            for k, ind in enumerate(rxn[species]):
                mol_reactant = self.mol_entries[ind]
                networkx_graph = mol_reactant.graph
                original_bonds_tmp = []
                # This is needed because there is a case where num_reactants != len(rxn['reactants']) or num_products != len(rxn['products'])
                if len(temp_species) <= k: 
                    break
                species_entry_ids.append(self.mol_entries[ind].entry_id)
                for i, j, weight in networkx_graph.edges:
                    original_bonds_tmp.append(sorted([i, j]))
                    species_total_bonds.add(tuple(sorted([temp_species[k][i], temp_species[k][j]])))
                original_bonds.append(original_bonds_tmp)
            return species_total_bonds, species_entry_ids, original_bonds
        
        reactants_total_bonds, reactants_entry_ids, original_bonds_react = find_total_bonds(rxn, 'reactants', reactants, products)
        products_total_bonds, products_entry_ids, original_bonds_prod = find_total_bonds(rxn, 'products', reactants, products)
        

        print(f"reactants_total_bonds: {reactants_total_bonds}")
        print(f"products_total_bonds: {products_total_bonds}")
        print(f"original bonds react: {original_bonds_react}")
        print(f"original bonds prod: {original_bonds_prod}")
        print(f"transformed atom map: {transformed_atom_map}")
        # find an union of bonds of reactants and products
        set_total_bonds = reactants_total_bonds.union(products_total_bonds)
        

        # convert to the correct format in "total_bonds" in mappings
        total_bonds = [[i,j] for i, j in set_total_bonds]

        # a dictionary, total_bonds_map, is used for creating bond_map for reactants and products in step 3
        total_bonds_map = {}
        for ind, bonds in enumerate(total_bonds):
            i, j = bonds
            total_bonds_map[(i,j)] = ind
        print(f"total_bonds_map: {total_bonds_map}")
        if rxn['is_redox']:
            assert len(set(reactants_total_bonds)) == len(set(products_total_bonds))
        
        # print(f'total_bonds: {total_bonds}')
        # print(f'atom_map: {atom_map}')
        # print(f'transformed_atom_map: {transformed_atom_map}')
        
            
        # step 3: Get bond_mapping
        bond_mapping = []

        # bond mapping for reactants
        bonds_in_reactants = [{} for _ in range(num_reactants)]
        for k, ind in enumerate(rxn['reactants']):
            mol_reactant = self.mol_entries[ind]
            networkx_graph = mol_reactant.graph
            # This is needed because there is a case where num_reactants != len(rxn['reactants']) or num_products != len(rxn['products'])
            if len(reactants) <= k: 
                break
            
            for bond_ind, edges in enumerate(networkx_graph.edges):
                i, j, _ = edges
                bonds_in_reactants[k][bond_ind] = total_bonds_map[tuple(sorted([reactants[k][i],reactants[k][j]]))]
        bond_mapping.append(bonds_in_reactants)

        # bond mapping for products
        bonds_in_products = [{} for _ in range(num_products)]
        for k, ind in enumerate(rxn['products']):
            mol_reactant = self.mol_entries[ind]
            networkx_graph = mol_reactant.graph
            # This is needed because there is a case where num_reactants != len(rxn['reactants']) or num_products != len(rxn['products'])
            if len(products) <= k: 
                break 
            for bond_ind, edges in enumerate(networkx_graph.edges):
                i, j, _ = edges
                bonds_in_products[k][bond_ind] = total_bonds_map[tuple(sorted([products[k][i],products[k][j]]))]
        bond_mapping.append(bonds_in_products)

        # print(f'bonds_in_reactants: {bonds_in_reactants}')
        # print(f'bonds_in_products: {bonds_in_products}')
        # print(f'reactants: {reactants}')
        # print(f'products: {products}')
        assert len(bonds_in_reactants) == len(reactants)
        assert len(bonds_in_products) == len(products)
        # for i in range(original_bonds_react):
        #     assert len(bonds_in_reactants[i]) == len(original_bonds_react[i])
        # for i in range(original_bonds_prod):
        #     assert len(bonds_in_products[i]) == len(original_bonds_prod[i])
        

        # step 4: get mapping
        mappings = {}
        mappings['bond_map'] = bond_mapping
        mappings['atom_map'] = transformed_atom_map
        mappings['total_bonds'] = total_bonds
        mappings['total_atoms'] = total_atoms
        mappings['num_bonds_total'] = len(total_bonds_map)
        mappings['num_atoms_total'] = len(total_atoms)
        if not rxn['is_redox']:
            print(f"mappings: {mappings}")
        #print(f"mapping: {mappings}")
        # print(f"atom_map: {atom_map}")
        # print(f"reactants: {reactants}")
        # print(f"products: {products}")

        # grab dgl graphs of reactants or products
        reactants_dgl_graphs  = [self.dgl_mol_dict[entry_i] for entry_i in reactants_entry_ids]       
        products_dgl_graphs = [self.dgl_mol_dict[entry_i] for entry_i in products_entry_ids]
        # print(f"reactants_dgl_graphs: {reactants_dgl_graphs}")
        # print(f"products_dgl_graphs: {products_dgl_graphs}")

        # create has_bonds in mappings
        has_bonds = defaultdict(list)
        for _ in range(len(reactants)):
            has_bonds['reactants'].append(True)
        for _ in range(len(products)):
            has_bonds['products'].append(True)

        # print(f"has_bonds: {has_bonds}")
        # print(f"mappings: {mappings}")

        # step 5: Create a reaction graphs and features
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

        # print(f"rxn_graph: {rxn_graph}")
        print(f"features: {features}")

        # step 5: update reaction features to the reaction graph
        for nt, ft in features.items():
            print(f"nt: {nt}")
            print(f"ft: {ft}")
            rxn_graph.nodes[nt].data.update({'ft': ft})

        # step 6: save a reaction graph and dG
        self.data[rxn_id] = {} # {'id': {}}
        self.data[rxn_id]['rxn_graph'] = rxn_graph
        self.data[rxn_id]['value'] = rxn['dG']  #torch.tensor([rxn['dG']])
        self.data[rxn_id]['reaction_features'] = features
        

        #### Write LMDB ####
        #1 load lmdb
        
        lmdb_path = "training_trial5.lmdb"
        lmdb_file = Path(lmdb_path)
        if lmdb_file.is_file():
            # file exists
            current_lmdb = LmdbDataset({'src': lmdb_path})
            lmdb_update = {
            "mean" : current_lmdb.mean,
            "std":   current_lmdb.std,
            "feature_size": self.grapher_features['feature_size'],
            "feature_name": self.grapher_features['feature_name'],
            "dtype": "float32"
            }
            current_length = current_lmdb.num_samples
        else:
            current_lmdb = lmdb.open(
                                    lmdb_path,
                                    map_size=1099511627776 * 2,
                                    subdir=False,
                                    meminit=False,
                                    map_async=True,
            
                                    )
            #2 initialize dict
            lmdb_update = {
            "mean" : 0,
            "std":   0,
            "feature_size": self.grapher_features['feature_size'],
            "feature_name": self.grapher_features['feature_name'],
            "dtype": "float32"
            }
            current_length = 0

     
        #3 update mean, std, feature_size, feature_name to dict_update
        #3.1 update mean
        prev_mean = lmdb_update["mean"]
        #print(f"current_length: {current_length}")
        n = current_length + 1
        current_y = rxn['dG']
        updated_mean = (current_y + (n-1)*prev_mean)/n
        lmdb_update["mean"] = updated_mean

        #3.2 update std
        prev_variance = (lmdb_update["std"])**2
        updated_variance = (n-1)/(n)*prev_variance + (n-1)/n*(prev_mean-updated_mean)**2 + (current_y - updated_mean)**2/n
        lmdb_update["std"] = math.sqrt(updated_variance)

        
        labels = {'value': torch.tensor([rxn['dG']]), 'value_rev': torch.tensor([0]), 'id': [str(rxn_id)], "reaction_type": ['']}
        data = (self.data[rxn_id]['rxn_graph'], self.data[rxn_id]['reaction_features'], labels)
        print(f"data: {data}")
        print(f"lmdb_update: {lmdb_update}")
        write_to_lmdb([data], current_length, lmdb_update, lmdb_path)

        

    def write_data(self): 
        # write a json file
        dumpfn(self.data, self.report_file_path)