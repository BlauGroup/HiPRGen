from HiPRGen.network_loader import *
from HiPRGen.species_questions import *
from HiPRGen.species_filter import *

network_loader = NetworkLoader(
    './scratch/li_test/rn.sqlite',
    './scratch/li_test/mol_entries.pickle',
    './scratch/li_test/initial_state.sqlite',
)

for m in network_loader.mol_entries:
    if m.has_covalent_ring:
        print(m)
