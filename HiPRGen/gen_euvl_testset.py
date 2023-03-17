log_message("Importing required files...")
from monty.serialization import loadfn, dumpfn
from mol_entry import MoleculeEntry
import initial_state
log_message("Done!")

log_message("Opening and reading species id file...")
with open('test_species.txt') as mol_id_doc: #open the text file and copy each line as the element of a list, removing '\n'
	mol_id_list = []
	for line in mol_id_doc:
		id_index = line.find("-")
		mol_id_list.append(line[id_index+1:len(line)-1])
log_message("Done!")

log_message("Opening json file...")
mol_json = "euvl_feb_full.json" #need to change to actually json
database_entries = loadfn(mol_json)
log_message("Done!")

log_message("Reading molecule entries from json file...")
if "has_props" in database_entries[0].keys():
        mol_entries= [MoleculeEntry.from_mp_doc(e) for e in database_entries]
else:
    mol_entries= [MoleculeEntry.from_dataset_entry(e) for e in database_entries]
log_message("Done!")

for molecule in mol_entries:
	target_mol_graph = MoleculeGraph.with_local_env_strategy(
        Molecule.from_file(xyz_file_path), OpenBabelNN()
    )

    # correction to the molecule graph
    target_mol_graph = metal_edge_extender(target_mol_graph)

    match = False
    index = -1
    while not match:
        index += 1
        mol_entry = mol_entries[index]
        species_mol_graph = mol_entry.mol_graph

        if mol_entry.charge == charge:
            match = target_mol_graph.isomorphic_to(species_mol_graph)

    if match:
        return mol_entry.ind
    else:
        return None

new_mol_entries = []

for elem in mol_id_list:
	ind = find_mol_entry_by_entry_id(mol_entries, elem)
	mol_entry = mol_entries[ind]
	mol_graph = mol_entry.mol_graph
	



