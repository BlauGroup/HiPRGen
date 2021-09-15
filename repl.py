from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *

network_loader = NetworkLoader(
    './scratch/li_test/rn.sqlite',
    './scratch/li_test/mol_entries.pickle',
    './scratch/li_test/initial_state.sqlite'
    )

reaction = {'number_of_reactants': 1, 'number_of_products': 1, 'reactants': (185, -1), 'products': (184, -1)}
