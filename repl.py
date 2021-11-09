from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *
from HiPRGen.network_renderer import *
import math

network_loader = NetworkLoader(
    '/home/danielbarter/HiPRGen/scratch/li_test/rn.sqlite',
    '/home/danielbarter/HiPRGen/scratch/li_test/mol_entries.pickle',
    '/home/danielbarter/HiPRGen/scratch/li_test/initial_state.sqlite',
    )


network_loader.load_trajectories()
network_loader.load_initial_state()
simulation_replayer = SimulationReplayer(network_loader)
pathfinder = Pathfinding(network_loader)

reactions_of_interest = {}
for sink in simulation_replayer.sinks:
    pathways = pathfinder.compute_pathways(sink)

    min_pathway = None
    weight = 1000000
    for pathway in pathways:
        if pathways[pathway]['weight'] < weight:
            weight = pathways[pathway]['weight']
            min_pathway = pathway

    for reaction_id in min_pathway:
        reactions_of_interest[reaction_id] = (184.0/256.0, 22/256.0, 22/256.0)


species_of_interest = {}
count = 0
for i in range(network_loader.number_of_species):
    if network_loader.initial_state_array[i] > 0:
        species_of_interest[i] = (0.5 + 0.47 * math.cos(math.pi + count * 0.01),
                                  0.5 + 0.47 * math.sin(math.pi + count * 0.01))
        count += 1

count = 0
for species_id in simulation_replayer.sinks:
    species_of_interest[species_id] = (0.5 + 0.47 * math.cos(count * 0.01),
                                       0.5 + 0.47 * math.sin(count * 0.01))
    count += 1



# network_loader = NetworkLoader(
#     '/home/danielbarter/big_network/rn.sqlite',
#     '/home/danielbarter/big_network/mol_entries.pickle',
#     '/home/danielbarter/big_network/initial_state.sqlite',
#     )




network_renderer = NetworkRenderer(
    network_loader,
    species_of_interest,
    reactions_of_interest,
    '/tmp/rn.png',
)

network_renderer.render()
