from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *
from HiPRGen.network_renderer import *

network_loader = NetworkLoader(
    '/home/danielbarter/HiPRGen/scratch/li_test/rn.sqlite',
    '/home/danielbarter/HiPRGen/scratch/li_test/mol_entries.pickle',
    '/home/danielbarter/HiPRGen/scratch/li_test/initial_state.sqlite',
    )


network_loader.load_trajectories()
network_loader.load_initial_state()
simulation_replayer = SimulationReplayer(network_loader)


reactions_of_interest = {}
for seed in network_loader.trajectories:
    for step in network_loader.trajectories[seed]:
        reaction_id = network_loader.trajectories[seed][step][0]
        if reaction_id not in reactions_of_interest:
            reactions_of_interest[reaction_id] = (184.0/256.0, 22/256.0, 22/256.0)


species_of_interest = {}
count = 0
for i in range(network_loader.number_of_species):
    if network_loader.initial_state_array[i] > 0:
        species_of_interest[i] = (0.01, 0.45 + count * 0.1)
        count += 1

count = 0
for species_id in simulation_replayer.sinks:
    species_of_interest[species_id] = (0.99, 0.1 + count * 0.1)
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
