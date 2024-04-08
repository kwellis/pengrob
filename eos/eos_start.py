"""
Equation of State Starters

A place to store functions that are used as starting points for various equation of states.
These functions provide a rough initial guess to start off the equation of state equations
"""

import math

import numpy as np


def whit_pk(mw_c7_plus: float) -> float:
    """Whitson Torp Pk

    Args:
        mw_c7_plus (float): Molecular Weight C7+, lb/lb-mol

    Returns:
        pk (float): Whitson Torp Pk
    """
    pk = 60 * mw_c7_plus - 4200
    return pk


def whit_capa(pabs: float, pk: float) -> float:
    """Whitson Torp A Factor

    Args:
        pabs (float): Evaluation Pressure, psia
        pk (float): Whitson Torp Pk, psia

    Returns:
        capa (float): Whitson Torp A Factor
    """
    capa = 1 - ((pabs - 14.7) / (pk - 14.7)) ** 0.6
    return capa


def whit_ki(pabs: float, tabs: float, pcrit: float, tcrit: float, acc: float, pk: float) -> float:
    """Whitson Torp Ki Factor

    Args:
        pabs (float): Evaluation Pressure, psia
        tabs (float): Evaluation Absolute Temperature, rankine
        pcrit (float): Critical Pressure, psia
        tcrit (float): Critical Temperature, rankine
        acc (float): Accentric Factor, unitless
        pk (float): Whitson Torp Pk, psia

    Returns:
        ki (float): Whitson Torp Ki
    """
    A = whit_capa(pabs, pk)
    ki = math.exp(5.37 * A * (1 + acc) * (1 - (tcrit / tabs))) * (pcrit / pabs) * (pcrit / pk) ** (A - 1)
    return ki


def wilson_ki(pabs: float, tabs: float, pcrit: float, tcrit: float, acc: float) -> float:
    """Wilson Ki Factor

    Args:
        pabs (float): Evaluation Absolute Pressure, psia
        tabs (float): Evaluation Absolute Temperature, rankine
        pcrit (float): Critical Pressure, psia
        tcrit (float): Critical Temperature, rankine
        acc (float): Accentric Factor, unitless

    Returns:
        ki (float): Wilson Ki
    """
    ki = math.exp(math.log(pcrit / pabs) + 5.373 * (1 + acc) * (1 - (tcrit / tabs)))
    return ki


def wilson_ki_list(pabs: float, tabs: float, ci_list: list, zi_list: list, prop_dict: dict, bini_dict: dict) -> list:
    """Wilson Equilibrium Constants

    Calculate the Wilson equilibrium constants for a mixture.

    Args:
        pabs (float): Absolute Evaluation Pressure, psia
        tabs (float): Absolute Evaluation Temp, rankine
        ci_list (list): List of String Components
        xi_list (list): Feed Phase Molar Fractions
        prop_dict (dict): Properties Dictionary
        bini_dict (dict): Binary Interaction Parameter Dictionary

    Returns:
        ki_list (float): Peng Robinson Equilibrium Constants
    """
    ki_list = []
    for ci, zi in zip(ci_list, zi_list):  # wilson is there to help with initial xi / yi split
        ki = wilson_ki(pabs, tabs, prop_dict[ci].pcrit, prop_dict[ci].tcrit, prop_dict[ci].acent)
        ki_list.append(ki)
    return ki_list


def safran_cfifteen(tabs: float, zi: float, pcrit: float, tcrit: float, acc: float) -> float:
    """Equation C-15 from Al-Safran Multiphase Flow

    Args:
        tabs (float): Evaluation Temperature Absolute, rankine
        zi (float): Molar Fraction of Component
        pcrit (float): Critical Pressure, psia
        tcrit (float): Critical Temperature, rankine
        acc (float): Accentric Factor, unitless

    Returns:
        pbi (float): Bubble Point Pressure for Component, psia
    """
    pbi = zi * pcrit * np.exp(5.371 * (1 + acc) * (1 - (tcrit / tabs)))
    return pbi


def bubblepoint_guess(tabs: float, comp_dict: dict, prop_dict: dict) -> float:
    """Bubble Point Guess

    Create a guessed value for bubblepoint, uses equation C-15 from Appendix of Multiphase
    Flow in piping by Al-Safran and Brill.

    Args:
        tabs (float): Evaluation Temperature Absolute, rankine
        comp_dict (dict): Mixture Molar Composition
        prop_dict (dict): Property Table for Lookup

    Returns:
        pbub (float): Guessed Bubblepoint Pressure, psia
    """
    plist = []

    for ci, zi in comp_dict.items():
        pbi = safran_cfifteen(tabs, zi, prop_dict[ci].pcrit, prop_dict[ci].tcrit, prop_dict[ci].acent)
        plist.append(pbi)
    return sum(plist)


def safran_ceighteen(tabs: float, zi: float, pcrit: float, tcrit: float, acc: float) -> float:
    """Equation C-18 from Al-Safran Multiphase Flow

    Args:
        tabs (float): Evaluation Temperature Absolute, rankine
        zi (float): Molar Fraction of Component
        pcrit (float): Critical Pressure, psia
        tcrit (float): Critical Temperature, rankine
        acc (float): Accentric Factor, unitless

    Returns:
        pdi (float): Dew Point Pressure for Component, psia
    """
    pdi = zi / (pcrit * np.exp(5.371 * (1 + acc) * (1 - (tcrit / tabs))))
    return pdi


def dewpoint_guess(tabs: float, comp_dict: dict, prop_dict: dict) -> float:
    """Dewpoint Guess

    Create a guessed value for dewpoint pressure, uses equation C-18 from Appendix of Multiphase
    Flow in piping by Al-Safran and Brill.

    Args:
        tabs (float): Evaluation Temperature Absolute, rankine
        comp_dict (dict): Mixture Molar Composition
        prop_dict (dict): Property Table for Lookup

    Returns:
        pdew (float): Guessed Dewpoint Pressure, psia
    """
    plist = []

    for ci, zi in comp_dict.items():
        pdi = safran_ceighteen(tabs, zi, prop_dict[ci].pcrit, prop_dict[ci].tcrit, prop_dict[ci].acent)
        plist.append(pdi)
    return 1 / sum(plist)
