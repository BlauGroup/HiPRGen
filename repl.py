from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *
from HiPRGen.network_renderer import *

# network_loader = NetworkLoader(
#     '/home/danielbarter/big_network/rn.sqlite',
#     '/home/danielbarter/big_network/mol_entries.pickle',
#     )

network_loader = NetworkLoader(
    '/home/danielbarter/HiPRGen/scratch/li_test/rn.sqlite',
    '/home/danielbarter/HiPRGen/scratch/li_test/mol_entries.pickle',
    )



network_renderer = NetworkRenderer(
    network_loader,
    None,
    None,
    '/tmp/rn.png',
)

network_renderer.render()
