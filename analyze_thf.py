import os
import sqlite3
import pickle
from HiPRGen.mc_analysis import *
from HiPRGen.network_loader import *
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge

base_dir = "/Users/ewcss/data/mg_network/thf_bcho_v4/high_potential"

network_loader = NetworkLoader(
    base_dir + '/rn.sqlite',
    base_dir + '/mol_entries.pickle'
    )

report_generator = ReportGenerator(
    network_loader.mol_entries,
    base_dir + '/main_report.tex',
    mol_pictures_folder_name=base_dir + '/mol_pictures',
    rebuild_mol_pictures=False)

h2_id = find_mol_entry_from_xyz_and_charge(
    network_loader.mol_entries,
    base_dir + "/mol_files/" +  'h2.xyz',
    0)

co_id = find_mol_entry_from_xyz_and_charge(
    network_loader.mol_entries,
    base_dir + "/mol_files/" +  'co.xyz',
    0)

c2h4_id = find_mol_entry_from_xyz_and_charge(
    network_loader.mol_entries,
    base_dir + "/mol_files/" +  'c2h4.xyz',
    0)

reaction_tally_report(
        network_loader,
        base_dir + '/reaction_tally.tex'
    )

pathfinding = Pathfinding(network_loader)
pathfinding.generate_pathway_report(
    h2_id,
    base_dir + '/h2_pathways.tex',
    sort_by_frequency=False
)

pathfinding.generate_pathway_report(
    co_id,
    base_dir + '/co_pathways.tex',
    sort_by_frequency=False
)

pathfinding.generate_pathway_report(
    c2h4_id,
    base_dir + '/c2h4_pathways.tex',
    sort_by_frequency=False
)
