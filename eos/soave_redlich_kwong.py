"""Soave-Redlich-Kwong Equation of State
"""

import numpy as np


def srk_mi(acc: float) -> float:
    """Soave-Redlich-Kwong m Factor

    Args:
        acc (float): Component Accentric Factor

    Returns:
        mi (float): SRK mi
    """
    mi = 0.48 + 1.574 * acc - 0.176 * acc**2
    return mi


def srk_alpha(tabs: float, tcrit: float, mi: float) -> float:
    """Soave-Redlich-Kwong Alpha

    Args:
        tabs (float): Evaluation Absolute Temperature, rankine
        tcrit (float): Critical Temperature, rankine
        mi (float): SRK mi

    Returns:
        alpha (float): SRK Alpha Factor
    """
    alpha = (1 + mi * (1 - (tabs / tcrit) ** 0.5)) ** 2
    return alpha


def srk_ai(pcrit: float, tcrit: float, rcon: float, alpha: float) -> float:
    """Soave-Redlich-Kwong ac Factor

    Args:
        pcrit (float): Critical Pressure, psia
        tcrit (float): Critical Temperature, rankine
        rcon (float): Gas Constant R, 10.731 blah blah
        alpha (float): SRK Alpha Factor

    Returns:
        ai (float): SRK ai, psia-ft6/(lbmol2)
    """
    ai = 0.42747 * alpha * rcon**2 * tcrit**2 / pcrit
    return ai


def srk_bi(pcrit: float, tcrit: float, rcon: float) -> float:
    """Soave-Redlich-Kwong b Factor

    Args:
        pcrit (float): Critical Pressure, psia
        tcrit (float): Critical Temperature, rankine
        rcon (float): Gas Constant R, 10.731 blah blah

    Returns:
        bi (float): SRK b
    """
    bi = 0.08664 * rcon * tcrit / pcrit
    return bi


def srk_capai(pabs: float, tabs: float, rcon: float, ai: float) -> float:
    """Soave-Redlich-Kwong A

    Args:
        pabs (float): Evaluation Absolute Pressure, psia
        tabs (float): Evaluation Absolute Temperature, rankine
        rcon (float): Gas Constant R, 10.731 blah blah
        ai (float): SRK ai, psia-ft6/(lbmol2)

    Returns:
        capa (float): SRK Capital A, unitless
    """
    capa = ai * pabs / (rcon**2 * tabs**2)
    return capa


def srk_capbi(pabs: float, tabs: float, rcon: float, bi: float) -> float:
    """Soave-Redlich-Kwong B

    Args:
        pabs (float): Evaluation Absolute Pressure, psia
        tabs (float): Evaluation Absolute Temperature, rankine
        rcon (float): Gas Constant R, 10.731 blah blah
        srk_b (float): SRK b Factor

    Returns:
        capb (float): SRK Capital B
    """
    capb = bi * pabs / (rcon * tabs)
    return capb


def srk_zfactors(A: float, B: float) -> np.ndarray:
    """Soave-Redlich-Kwong Z Factors

    Find the 3 roots of the Soave-Redlich-Kwong Equation of state.
    Return only the list of real roots, discard imaginary.

    Args:
        A (float): Soave-Redlich-Kwong A
        B (float): Soave-Redlich-Kwong B

    Return:
        zray (np array): List of real Z Factors
    """
    coeff = [1, -1, A - B + B**2, -A * B]  # coefficients of srk cubic equation
    zroots = np.roots(coeff)
    zray = zroots.real[abs(zroots.imag) < 1e-5]  # filter out imaginary numbers
    return zray
