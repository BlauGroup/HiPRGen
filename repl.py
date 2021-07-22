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
    'reactants' : [29, 161],
    'products' : [32, 160]
}

fragment_matching_found(reaction, network_loader.mol_entries, params)
