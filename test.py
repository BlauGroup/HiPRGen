import os
import sqlite3
import pickle
from HiPRGen.mc_analysis import *
from HiPRGen.network_loader import *
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge


class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'

network_loader = NetworkLoader(
    './scratch/rn.sqlite',
    './scratch/ronald_mol_entries.pickle'
    )


def test_species_filter(network_loader):
    with open('./data/ronald_mol_entries.pickle', 'rb') as f:
        mol_entries_data = pickle.load(f)



    if network_loader.mol_entries == mol_entries_data:
        print(bcolors.PASS + "passed: test_species_filter" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "failed: test_species_filter" + bcolors.ENDC)

def test_reaction_gen(network_loader):

    if network_loader.number_of_reactions == 5525:
        print(bcolors.PASS + "passed: test_reaction_gen" + bcolors.ENDC)
    else:
        print(bcolors.FAIL + "failed: test_reaction_gen" + bcolors.ENDC)


def test_mc_pathfinding(network_loader):

    report_generator = ReportGenerator(
        network_loader.mol_entries,
        './scratch/dummy.tex',
        rebuild_mol_pictures=True)

    reaction_tally_report(
        network_loader,
        './scratch/reaction_tally.tex'
    )

    LEDC_id = find_mol_entry_from_xyz_and_charge(
        network_loader.mol_entries,
        './xyz_files/LEDC.xyz',
        0)

    pathfinding = Pathfinding(network_loader)
    pathfinding.generate_pathway_report(
        LEDC_id,
        './scratch/LEDC_pathways.tex',
        sort_by_frequency=False
    )

    species_report(network_loader, './scratch/species_report.tex')
    sink_report(network_loader, './scratch/sink_report.tex')


test_mc_pathfinding(network_loader)
test_species_filter(network_loader)
test_reaction_gen(network_loader)
