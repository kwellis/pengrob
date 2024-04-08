"""Numerical Methods

A place to keep various optimization and numerical iteration routies.
"""

import math


# https://docs.python.org/3/library/string.html#formatspec
def mix_comp_table(mix_comp: dict, xi_list: list, yi_list: list) -> str:
    """Creates a fancy table to see some of the stored data"""
    sformat = "{:>8} | {:>8} | {:>8} | {:>8} \n"  # string format
    nformat = "{:>8} | {:>8.4f} | {:>8.4f} | {:>8.4f} \n"  # number format
    spc = 48 * "-" + "\n"  # spacing
    pout = sformat.format("Hydro", "Overall", "Liquid", "Vapor")
    pout = pout + sformat.format("Carbon", "mol frac", "mol frac", "mol frac") + spc
    for name, zi, xi, yi in zip(mix_comp.keys(), mix_comp.values(), xi_list, yi_list):
        pout += nformat.format(name, zi, xi, yi)
    return pout


def psi_secant(psi1: float, psi2: float, ysum1: float, ysum2: float) -> float:
    """Next Bubble Point Pressure with Secant Method

    Uses the secant method to calculate the next psi near bubblepoint.

    Args:
        psi1 (float): Pressure One, psig
        psi2 (float): Pressure Two, psig
        ysum1 (float): Liquid Phase Summation at Psi1
        ysum2 (float): Liquid Phase Summation at Psi2

    Return:
        psi3 (float): Pressure Three, psig
    """
    ysum1 -= 1  # subtract one, to solve to zero
    ysum2 -= 1  # subtract one, to solve to zero
    psi3 = psi2 - ysum2 * (psi1 - psi2) / (ysum1 - ysum2)
    return psi3


def nv_secant(nv1: float, nv2: float, xsum1: float, xsum2: float) -> float:
    """Next Bubble Point Pressure with Secant Method

    Uses the secant method to calculate the next nv alue.

    Args:
        nv1 (float): Vapor Mole Fraction One, psig
        nv2 (float): Vapor Mole Fraction Two, psig
        xsum1 (float): Liquid Phase Summation at Psi1
        xsum2 (float): Liquid Phase Summation at Psi2

    Return:
        nv3 (float): Pressure Three, psig
    """
    xsum1 -= 1  # subtract one, to solve to zero
    xsum2 -= 1  # subtract one, to solve to zero
    nv3 = nv2 - xsum2 * (nv1 - nv2) / (xsum1 - xsum2)
    return nv3
