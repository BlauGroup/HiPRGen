from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *

network_loader = NetworkLoader(
    './scratch/rn.sqlite',
    './scratch/ronald_mol_entries.pickle'
    )

params={
    'temperature' : ROOM_TEMP,
    'electron_free_energy' : -1.4
    }

reaction = {
    'number_of_reactants' : 1,
    'number_of_products' : 1,
    'reactants' : [54],
    'products' : [53]
}

m54 = network_loader.mol_entries[54]
m53 = network_loader.mol_entries[53]
