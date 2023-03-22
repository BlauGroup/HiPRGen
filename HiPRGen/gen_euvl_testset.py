print("Importing required files...")
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from monty.serialization import loadfn, dumpfn
from mol_entry import MoleculeEntry
from initial_state import find_mol_entry_by_entry_id
import copy
print("Done!")

print("Opening and reading species id file...")
with open('test_species.txt') as mol_id_doc: #open the text file and copy each line as the element of a list, removing '\n'
    mol_id_list = []
    for line in mol_id_doc:
        if line != '':
            id_index = line.find("-")
            mol_id_list.append(line[id_index+1:len(line)-1])
print("Done!")

print("Opening json file...")
mol_json = "euvl_feb_full.json"
database_entries = loadfn(mol_json) #loads the json file to a dictionary
print("Done!")

print("Reading molecule entries from json file...") #converts dictionary from loaded json file to a list
mol_entries= [MoleculeEntry.from_mp_doc(e) for e in database_entries]
print("Done!")


final_mol_ids = []

print('copying desired molecules to a new list...')
for mol_id in mol_id_list: #takes the list of ids, finds their corresponding molecules in the json file (as well as differently charged variants) and adds them to the list
    for m in mol_entries:
        m_id = m.entry_id
        if m_id == mol_id:
            graph = m.mol_graph
            for mol_entry in mol_entries:
                molecule_graph = mol_entry.mol_graph #test if other charges found
                if molecule_graph.isomorphic_to(graph):
                    final_mol_ids.append(mol_entry.entry_id)
print('Done!')
print(len(final_mol_ids)) #currently missing 3 entries

test_set = []
for e in database_entries:
    if e['molecule_id'] in final_mol_ids:
        entry = copy.deepcopy(e)
        val = entry.pop('nbo_population')
        val = entry.pop('nbo_lone_pairs')
        val = entry.pop('nbo_bonds')
        val = entry.pop('nbo_interactions')
        val = entry.pop('alpha_population')
        val = entry.pop('beta_population')
        val = entry.pop('alpha_lone_pairs')
        val = entry.pop('beta_lone_pairs')
        val = entry.pop('alpha_bonds')
        val = entry.pop('beta_bonds')
        val = entry.pop('alpha_interactions')
        val = entry.pop('beta_interactions')
        test_set.append(entry)
        
print(len(test_set))
dumpfn(test_set, 'euvl_test_set.json')