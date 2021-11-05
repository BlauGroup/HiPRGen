from HiPRGen.network_loader import *
from HiPRGen.reaction_questions import *
from HiPRGen.mc_analysis import *
from matplotlib.colors import hsv_to_rgb

network_loader = NetworkLoader(
    './scratch/li_test/rn.sqlite',
    './scratch/li_test/mol_entries.pickle',
    './scratch/li_test/initial_state.sqlite'
    )



simulation_replayer = SimulationReplayer(network_loader)
simulation_replayer.time_series_graph(
    network_loader.trajectories.keys(),
    simulation_replayer.sinks,
    '/tmp/avg.svg'
)
