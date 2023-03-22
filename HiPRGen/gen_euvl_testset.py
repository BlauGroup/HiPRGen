print("Importing required files...")
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from monty.serialization import loadfn, dumpfn
from mol_entry import MoleculeEntry
from initial_state import find_mol_entry_by_entry_id
print("Done!")

print("Opening and reading species id file...")
with open('test_species.txt') as mol_id_doc: #open the text file and copy each line as the element of a list, removing '\n'
    mol_id_list = []
    for line in mol_id_doc:
        id_index = line.find("-")
        mol_id_list.append(line[id_index+1:len(line)-1])
print("Done!")

print("Opening json file...")
mol_json = "euvl_feb_full.json"
database_entries = loadfn(mol_json)
print("Done!")

print("Reading molecule entries from json file...")
if "has_props" in database_entries[0].keys():
        mol_entries= [MoleculeEntry.from_mp_doc(e) for e in database_entries]
else:
    mol_entries= [MoleculeEntry.from_dataset_entry(e) for e in database_entries]
print("Done!")


new_mol_entries = []

print('copying desired molecules to a new list...')
for mol_id in mol_id_list: #takes the list of ids, finds their corresponding molecules in the json file (as well as differently charged variants) and adds them to the list
    for m in mol_entries:
        m_id = m.entry_id
        if m_id == mol_id:
            graph = m.mol_graph
            for molecule in mol_entries:
                molecule_graph = molecule.mol_graph
                if molecule_graph.isomorphic_to(graph):
                    new_mol_entries.append(molecule)
print('Done!')
print(new_mol_entries)

mol_entries_dict = {}
for i in new_mol_entries: # All molecules formed at least once are in sink_data
    mol_entries_dict[i.entry_id] = i
print('Dumping desired molecules to a new json file...')
dumpfn(new_mol_entries, 'euvl_test_set.json')