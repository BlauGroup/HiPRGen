from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *
from HiPRGen.network_renderer import *
import math

network_loader = NetworkLoader(
    './scratch/li_test/rn.sqlite',
    './scratch/li_test/mol_entries.pickle',
    './scratch/li_test/initial_state.sqlite',
    )

network_loader.load_trajectories()
network_loader.load_initial_state()
pathfinding = Pathfinding(network_loader)
simulation_replayer = SimulationReplayer(network_loader)

reactions_which_fired = set()
for seed in network_loader.trajectories:
    for step in network_loader.trajectories[seed]:
        reaction_id = network_loader.trajectories[seed][step][0]
        reactions_which_fired.add(reaction_id)



species_of_interest = {}
count = 0
for i in range(network_loader.number_of_species):
    if network_loader.initial_state_array[i] > 0:
        species_of_interest[i] = (0.5 + 0.47 * math.cos(math.pi + count * 0.02),
                                  0.5 + 0.47 * math.sin(math.pi + count * 0.02))
        count += 1



count = 0
for i in simulation_replayer.sinks:
    species_of_interest[i] = (0.5 + 0.47 * math.cos(count * 0.02),
                              0.5 + 0.47 * math.sin(count * 0.02))
    count += 1


top_pathway_reactions = set()
for reaction_id in simulation_replayer.sinks:
    pathways = pathfinding.compute_pathways(reaction_id)
    pathways_sorted = sorted(pathways, key=lambda p: pathways[p]['weight'])[0:4]
    for p in pathways_sorted:
        top_pathway_reactions.update(p)




network_renderer = NetworkRenderer(
    network_loader,
    species_of_interest,
    reactions_which_fired,
    '/tmp/rn_reactions_which_fired.png',
)

network_renderer = NetworkRenderer(
    network_loader,
    species_of_interest,
    top_pathway_reactions,
    '/tmp/top_pathway_reactions.png',
)

network_renderer.render()
