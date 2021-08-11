# coding: utf-8
# Copyright (c) MR.Net development team


from enum import Enum

# Basic constants

# Room temperature (25 C) in Kelvin
ROOM_TEMP = 298.15

# Boltzmann constant in eV / K
KB = 8.617333262 * 10 ** -5

# Planck constant in eV * s
PLANCK = 4.135667696 * 10 ** -15

class Terminal(Enum):
    KEEP = 1
    DISCARD = -1

metals = frozenset(["Li", "Na", "K", "Mg", "Ca", "Zn", "Al"])
m_formulas = frozenset([m + "1" for m in metals])


li_solvation_correction = {
    "Li" : -0.68,
}

g2_solvation_correction = {
    "Mg_1": -0.56,
    "Mg_2": -1.49
}

thf_solvation_correction = {
    "Mg_1": -0.70,
    "Mg_2": -1.91
}


li_coordination_radius = {
    "Li" : 2.4
}

mg_coordination_radius = {
    "Mg_1": 2.4,
    "Mg_2": 2.4
}


li_max_number_of_coordination_bonds = {
    "Li" : 4
}

mg_max_number_of_coordination_bonds = {
    "Mg_1": 5,
    "Mg_2": 6
}
