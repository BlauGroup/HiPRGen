import sys
from functools import partial
from monty.serialization import loadfn
from HiPRGen.species_filter import *
from HiPRGen.bucketing import *
from HiPRGen.report_generator import ReportGenerator
from HiPRGen.species_questions import *
from HiPRGen.solvation_corrections import (solvation_correction_ec,
                                           solvation_correction_g2,
                                           solvation_correction_thf)

if len(sys.argv) != 5:
    print("usage: python run_bucketing.py database_json_file mol_entries_pickle_file bucket_db_file report_file")
    quit()

database_json_file = sys.argv[1]
mol_entries_pickle_file = sys.argv[2]
bucket_db_file = sys.argv[3]
report_file = sys.argv[4]

decision_tree = [
    (fix_hydrogen_bonding_length, Terminal.KEEP),
    (metal_ion_filter, Terminal.DISCARD),
    (bad_metal_coordination, Terminal.DISCARD),
    (mol_not_connected, Terminal.DISCARD),
    (metal_complex, Terminal.DISCARD),
    (add_star_hashes, Terminal.KEEP),

    (partial(add_unbroken_fragment,
             fragment_neighborhood_width), Terminal.KEEP),
    (partial(add_single_bond_fragments,
             fragment_neighborhood_width), Terminal.KEEP),

    (default_true, Terminal.KEEP)
    ]

database_entries = loadfn(database_json_file)
mol_entries = species_filter(
    database_entries,
    mol_entries_pickle_file,
    report_file,
    solvation_correction=None,
    species_decision_tree=decision_tree,
    generate_unfiltered_mol_pictures=False
)
bucket(mol_entries, bucket_db_file)

