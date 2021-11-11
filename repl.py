from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *
from HiPRGen.network_renderer import *
import math
from multiprocessing import Pool

network_loader = NetworkLoader(
    './scratch/li_test/rn.sqlite',
    './scratch/li_test/mol_entries.pickle',
    './scratch/li_test/initial_state.sqlite',
    )

network_loader.load_trajectories()
network_loader.load_initial_state()
pathfinding = Pathfinding(network_loader)
simulation_replayer = SimulationReplayer(network_loader)

render_reactions_which_fired(network_loader, '/tmp/reactions_which_fired.png')
