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
    'number_of_reactants' : 2,
    'number_of_products' : 2,
    'reactants' : [53,53],
    'products' : [10,191]
}

m53 = network_loader.mol_entries[53]
m191 = network_loader.mol_entries[191]
f0=m53.fragment_data[4].fragments[1]
f1=m191.fragment_data[0].fragments[0]
