# coding: utf-8
# Copyright (c) MR.Net development team


from enum import Enum
from monty.json import MSONable

# Basic constants

# Room temperature (25 C) in Kelvin
ROOM_TEMP = 298.15

# Boltzmann constant in eV / K
KB = 8.617333262 * 10 ** -5

# Planck constant in eV * s
PLANCK = 4.135667696 * 10 ** -15

class Terminal(MSONable, Enum):
    KEEP = 1
    DISCARD = -1

metals = frozenset(["Li", "Na", "K", "Mg", "Ca", "Zn", "Al"])
m_formulas = frozenset([m + "1" for m in metals])


# solvation environments
li_ec = {
    "solvation_correction" : {
        "Li" : -0.68
    },

    "coordination_radius" : {
        "Li" : 2.4
    },

    "max_number_of_coordination_bonds" : {
        "Li" : 4
    }
}


mg_g2 = {
    "solvation_correction" : {
        "Mg_1": -0.56,
        "Mg_2": -1.49
    },

    "coordination_radius" : {
        "Mg_1": 2.4,
        "Mg_2": 2.4
    },

    "max_number_of_coordination_bonds" : {
        "Mg_1": 5,
        "Mg_2": 6
    }
}


mg_thf = {
    "solvation_correction" : {
        "Mg_1": -0.70,
        "Mg_2": -1.91
    },

    "coordination_radius" : {
        "Mg_1": 2.4,
        "Mg_2": 2.4
    },

    "max_number_of_coordination_bonds" : {
        "Mg_1": 5,
        "Mg_2": 6
    }
}
