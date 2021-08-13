from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *

network_loader = NetworkLoader(
    './scratch/mg_test/rn.sqlite',
    './scratch/mg_test/mol_entries.pickle'
    )

params={
    'temperature' : ROOM_TEMP,
    'electron_free_energy' : -2.06
    }

reaction = {
    'number_of_reactants' : 2,
    'number_of_products' : 1,
    'reactants' : [43,89],
    'products' : [36]
}

