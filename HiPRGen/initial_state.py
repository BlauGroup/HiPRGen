from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.fragmenter import metal_edge_extender
import sqlite3


def find_mol_entry_from_xyz_and_charge(mol_entries, xyz_file_path, charge):
    """
    given a file 'molecule.xyz', find the mol_entry corresponding to the
    molecule graph with given charge
    """
    target_mol_graph = MoleculeGraph.with_local_env_strategy(
        Molecule.from_file(xyz_file_path), OpenBabelNN()
    )

    # correction to the molecule graph
    target_mol_graph = metal_edge_extender(target_mol_graph)

    match = False
    index = -1
    while not match:
        index += 1
        mol_entry = mol_entries[index]
        species_mol_graph = mol_entry.mol_graph

        if mol_entry.charge == charge:
            match = target_mol_graph.isomorphic_to(species_mol_graph)

    if match:
        return mol_entry.ind
    else:
        return None

def find_mol_entry_by_entry_id(mol_entries, entry_id):
    """
    given an entry_id, return the corresponding mol enentry index
    """

    for m in mol_entries:
        if m.entry_id == entry_id:
            return m.ind

create_initial_state_table = """
    CREATE TABLE initial_state (
            species_id             INTEGER NOT NULL PRIMARY KEY,
            count                  INTEGER NOT NULL
    );
"""

create_trajectories_table = """
    CREATE TABLE trajectories (
            seed         INTEGER NOT NULL,
            step         INTEGER NOT NULL,
            reaction_id  INTEGER NOT NULL,
            time         REAL NOT NULL
    );
"""

create_factors_table = """
    CREATE TABLE factors (
            factor_zero      REAL NOT NULL,
            factor_two       REAL NOT NULL,
            factor_duplicate REAL NOT NULL
    );
"""


def insert_initial_state(
        initial_state,
        mol_entries,
        initial_state_db,
        factor_zero = 1.0,
        factor_two = 1.0,
        factor_duplicate = 0.5
):
    """
    initial state is a dict mapping species ids to counts.
    """

    rn_con = sqlite3.connect(initial_state_db)
    rn_cur = rn_con.cursor()
    rn_cur.execute(create_initial_state_table)
    rn_cur.execute(create_trajectories_table)
    rn_cur.execute(create_factors_table)
    rn_con.commit()

    rn_cur.execute(
        "INSERT INTO factors VALUES (?,?,?)",
        (factor_zero, factor_two, factor_duplicate))

    num_species = len(mol_entries)


    for i in range(num_species):
        rn_cur.execute(
            "INSERT INTO initial_state VALUES (?,?)",
            (i, initial_state.get(i,0)))

    rn_con.commit()



