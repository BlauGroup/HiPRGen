from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *
from HiPRGen.network_renderer import *
import sqlite3
import math

gold = (0, 0, 0)
green = (0/256, 103/256, 0/256)
light_green = (132/256, 176/256, 0/256)
purple = (149/256, 63/256, 202/256)
black = (0,0,0)

network_loader = NetworkLoader(
    '/global/scratch/dbarter/oct_production_networks/li_1.4/rn.sqlite',
    '/global/scratch/dbarter/oct_production_networks/li_1.4/mol_entries.pickle',
    '/global/scratch/dbarter/oct_production_networks/li_1.4/li_ec_co2/initial_state.sqlite',
)

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







# network_loader = NetworkLoader(
#     './scratch/li_test/rn.sqlite',
#     './scratch/li_test/mol_entries.pickle',
#     './scratch/li_test/initial_state.sqlite',
# )

# colors = {
#     75 : green,
#     88 : black,
#     58 : purple,
#     183 : green,
#     63 : green,
#     6 : black,
#     1 : green
# }

network_loader.load_trajectories()
network_loader.load_initial_state()




pathfinding = Pathfinding(network_loader)
simulation_replayer = SimulationReplayer(network_loader)

render_species(network_loader, '../network_renders/frame_1.png')


render_reactions_which_fired(
    network_loader,
    [],
    '../network_renders/frame_2.png'
)

render_reactions_which_fired(
    network_loader,
    colors,
    '../network_renders/frame_3.png'
)

render_reactions_which_fired_new_positions(
    network_loader,
    colors,
    '../network_renders/frame_4.png')

# clear out the database connections so that we can pass the pathfinding
# objects into Pool(n).map
network_loader.rn_con = None
network_loader.initial_state_con = None
network_loader.cur = None


render_top_pathways(
    pathfinding,
    colors,
    '../network_renders/frame_5.png',
    num_threads=20
)


render_top_highlighted(
    pathfinding,
    colors,
    '../network_renders/frame_6.png',
    58,
    num_threads=20
)


