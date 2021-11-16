from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *
from HiPRGen.network_renderer import *
import sqlite3
import math


network_loader = NetworkLoader(
    './scratch/li_test/rn.sqlite',
    './scratch/li_test/mol_entries.pickle',
    './scratch/li_test/initial_state.sqlite',
)


network_loader.load_trajectories()
network_loader.load_initial_state()


gold = (222/256, 146/256, 4/256)
green = (13/256, 79/256, 18/256)
light_green = (4/256, 222/256, 0/256)
purple = (86/256, 10/256, 110/256)
black = (0,0,0)

colors = {
    # Gases:
    4048 : green,
    1239 : green,

    # Inorganics:
    2539 : green,
    3418 : green,
    1459 : green,
    4432 : green,

    # Alkyl carbonates:
    4329 : green,
    23 : light_green,
    2549 : green,
    3867 : light_green,
    173 : green,
    496 : green,
    3827 : green,

    # Carbs, esters, oxides:
    2711 : black,
    2750 : black,

    # Cyclic species:
    2341 : gold,
    2676 : purple,
    3829 : black,
}

pathfinding = Pathfinding(network_loader)
simulation_replayer = SimulationReplayer(network_loader)


render_reactions_which_fired(
    network_loader,
    simulation_replayer.sinks,
    colors,
    '/tmp/reactions_which_fired.png'
)

network_loader.rn_con = None
network_loader.initial_state_con = None
network_loader.cur = None


render_top_pathways(
    pathfinding,
    simulation_replayer.sinks,
    colors,
    '/tmp/top_pathways.png'
)


