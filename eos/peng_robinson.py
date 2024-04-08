"""Peng Robinson Equation of State

https://www.e-education.psu.edu/png520/m16_p6.html
https://www.e-education.psu.edu/png520/m17.html
"""

import math

import numpy as np

import eos.mixing_rules as mr


def pengrob_mi(acc: float) -> float:
    """Peng Robinson mi Factor

    Args:
        acc (float): Accentric Factor, unitless

    Returns:
        mi (float): Peng Robinson mi
    """
    if acc < 0.49:
        mi = 0.37464 + 1.54226 * acc - 0.26922 * acc**2
    else:  # robinson correction
        mi = 0.3796 + 1.485 * acc - 0.1644 * acc**2 + 0.01667 * acc * 3
    return mi


def pengrob_alpha(tabs: float, tcrit: float, mi: float) -> float:
    """Peng Robinson Alpha

    Args:
        tabs (float): Eval Absolute Temperature, rankine
        tcrit (float): Critical Temperature, rankine
        mi (float): Peng Robinson mi

    Returns:
        alpha (float): Peng Robinson Alpha Value
    """
    alpha = (1 + mi * (1 - (tabs / tcrit) ** (1 / 2))) ** 2
    return alpha


def pengrob_ai(pcrit: float, tcrit: float, rcon: float, alpha: float) -> float:
    """Peng Robinson Component ai

    Args:
        pcrit (float): Critical Pressure, psia
        tcrit (float): Critical Temperature, rankine
        rcon (float): Gas Constant R, 10.731 psia-ft3/(lbmol-R)
        alpha (float): Peng Robinson Alpha Value

    Returns:
        ai (float): Peng Robinson ai, psia-ft6/(lbmol2)
    """
    ai = 0.45724 * alpha * rcon**2 * tcrit**2 / pcrit
    return ai


def pengrob_bi(pcrit: float, tcrit: float, rcon: float) -> float:
    """Peng Robinson Component bi

    Args:
        pcrit (float): Critical Pressure, psia
        tcrit (float): Critical Temperature, rankine
        rcon (float): Gas Constant R, 10.731 psia-ft3/(lbmol-R)

    Returns:
        bi (float): Peng Robinson bi, ft3/lbmol
    """
    bi = 0.0778 * rcon * tcrit / pcrit
    return bi


def pengrob_capai(pabs: float, tabs: float, rcon: float, ai: float) -> float:
    """Peng Robinson Capital Ai

    Args:
        pabs (float): Eval Absolute Pressure, psia
        tabs (float): Eval Absolute Temperature, rankine
        rcon (float): Gas Constant R, 10.731 psia-ft3/(lbmol-R)
        ai (float): Peng Robinson ai, psia-ft6/(lbmol2)

    Returns:
        capai (float): Peng Robinson Ai, unitless?
    """
    capai = ai * pabs / (rcon**2 * tabs**2)
    return capai


def pengrob_capbi(pabs: float, tabs: float, rcon: float, bi: float) -> float:
    """Peng Robinson Capital Bi

    Args:
        pabs (float): Eval Absolute Pressure, psia
        tabs (float): Eval Absolute Temperature, rankine
        rcon (float): Gas Constant R, 10.731 psia-ft3/(lbmol-R)
        bi (float): Peng Robinson bi, ft3/lbmol

    Returns:
        capbi (float): Peng Robinson Bi, unitless?
    """
    capbi = bi * pabs / (rcon * tabs)
    return capbi


# should I make this decide what to return for liquid, gas or single phase?
def pengrob_zfactors(A: float, B: float) -> np.ndarray:
    """Peng Robinson Z Factors

    Find the 3 roots of the Peng Robinson Equation of state.
    Return only the list of real roots, discard imaginary.

    Args:
        A (float): Peng Robinson A
        B (float): Peng Robinson B

    Return:
        zray (np array): List of real Z Factors
    """
    coeff = [1, -1 * (1 - B), A - 2 * B - 3 * B**2, -1 * (A * B - B**2 - B**3)]  # coefficients of cubic equation
    zroots = np.roots(coeff)
    zray = zroots.real[abs(zroots.imag) < 1e-5]  # filter out imaginary numbers
    return zray


def pengrob_ab_rays(tabs: float, comp_list: list, prop_dict: dict) -> tuple[list, list]:
    """Peng Robinson a and b Arrays

    Create two lists of Peng Robinson a and b values for components.
    All you actually need here is a list of components you are using.
    You don't actually need the molar fractions at this point.

    Args:
        tabs (float): Absolute Temperature, Rankine
        comp_list (list): Feed Components, string of values
        prop_dict (dict): Critical Property Lookup Dictionary

    Returns:
        pra_list (list): Peng Robinson a values for each component in feed mixture
        prb_list (list): Peng Robinson b values for each component in feed mixture
    """
    # should I just hard code rcon everywhere else? make a class?
    rcon = 10.731  # psia-ft3/(lbmol-R)
    pra_list = []
    prb_list = []
    # comp is component string
    for comp in comp_list:
        mi = pengrob_mi(prop_dict[comp].acent)
        pcrit = prop_dict[comp].pcrit  # reduce number of dictionary lookups
        tcrit = prop_dict[comp].tcrit
        alpha = pengrob_alpha(tabs, tcrit, mi)
        ai = pengrob_ai(pcrit, tcrit, rcon, alpha)
        bi = pengrob_bi(pcrit, tcrit, rcon)
        pra_list.append(ai)
        prb_list.append(bi)
    return pra_list, prb_list


def pengrob_fugj(ci: str, cj_list: list, zj_list: list, aj_list: list, bini_dict: dict) -> float:
    """Peng Robinson Fugacity Coefficient Summation

    Completes the summation that is inside the Peng Robinson Fugacity Coefficient
    for each component to be evaluated. I actually like the term fugj better, but
    I don't have a great way to rename stuff? I guess sum is not worst than list...

    Args:
        ci (str): Primary Component Abbreviation
        cj_list (list): Secondary Components in the Mixture
        zj_list (list): Molar Fractions of the Components
        aj_list (list): Peng Rob a values for Secondary Components
        bini_dict (dict): Binary Interaction Parameters for Lookup

    Returns:
        fugsum (float): Summation of j components in PR Fugacity Coefficient
    """
    fugj_list = []
    for cj, zj, aj in zip(cj_list, zj_list, aj_list):
        kij = bini_dict[ci][cj]
        fugj = zj * (1 - kij) * (aj) ** (1 / 2)
        fugj_list.append(fugj)
    fugsum = sum(fugj_list)
    return fugsum


def pengrob_fugco(
    ai: float,
    bi: float,
    am: float,
    bm: float,
    Am: float,
    Bm: float,
    Zm: float,
    fugj: float,
) -> float:
    """Peng Robinson Fugacity Coefficient

    Calculate the fugacity coefficient for a single component in either
    the vapor or liquid phase.

    Args:
        ai (float): Peng Robinson a for Component
        bi (float): Peng Robinson b for Component, ft3/lbmol
        am (float): Peng Robinson a for Mixture
        bm (float): Peng Robinson b for Mixture
        Am (float): Peng Robinson A for Mixture
        Bm (float): Peng Robinson B for Mixture
        Zm (float): Z Factor Overall for Vapor or Liquid
        fugj (float): j Component Summation from Peng Robinson Fug Coeff.

    Returns:
        phi_i (float): Peng Robinson Fugacity Coefficient for Comp i
    """
    # end of Peng Robinson Fugacity Equation
    fugend = math.log((Zm + (math.sqrt(2) + 1) * Bm) / (Zm - (math.sqrt(2) - 1) * Bm))

    # you get a super small imaginary value here...?
    lnphi = (
        -math.log(Zm - Bm)
        + (Zm - 1) * bi / bm  # noqa: W503
        - Am / (2**1.5 * Bm) * ((1 / am) * (2 * math.sqrt(ai) * fugj) - bi / bm) * fugend  # noqa: W503
    )
    phi_i = math.exp(lnphi)
    return phi_i


def pengrob_fugco_list(
    pabs: float, tabs: float, comp_list: list, zi_list: list, prop_dict: dict, bini_dict: dict, vapor: bool
) -> list:
    """Peng Robinson Fugacity Coefficient List for Liquids or Vapor

    Calculate the fugacity coefficients for liquid or Vapor phase.

    Args:
        pabs (float): Absolute Evaluation Pressure, psia
        tabs (float): Absolute Evaluation Temp, rankine
        comp_list (list): List of String Components
        zi_list (list): Molar Fractions of Evaluated Mixture
        prop_dict (dict): Properties Dictionary
        bini_dict (dict): Binary Interaction Parameter Dictionary
        vapor (bool): True - Evaluate Vapor, False Evaluate Liquid

    Returns:
        phi_list (float): Peng Robinson Fugacity Coefficient for specified phase
    """
    rcon = 10.731  # psia-ft3/lbmol-R
    ai_list, bi_list = pengrob_ab_rays(tabs, comp_list, prop_dict)  # same for both

    # make mixing rules a seperate file? also make one for loop for both?
    amix = mr.mix_a(comp_list, zi_list, ai_list, bini_dict)
    bmix = mr.mix_b(zi_list, bi_list)

    Amix = pengrob_capai(pabs, tabs, rcon, amix)
    Bmix = pengrob_capbi(pabs, tabs, rcon, bmix)

    zray = pengrob_zfactors(Amix, Bmix)

    if vapor:
        zfac = max(zray)  # keep the largest value for vapor
    else:
        zfac = min(zray)  # keep smallest value for liquid

    phi_list = []

    aj_list = ai_list.copy()  # prevent weird behavior from occuring with loop
    cj_list = comp_list.copy()

    for ci, ai, bi in zip(comp_list, ai_list, bi_list):
        fugj = pengrob_fugj(ci, cj_list, zi_list, aj_list, bini_dict)
        phi_i = pengrob_fugco(ai, bi, amix, bmix, Amix, Bmix, zfac, fugj)
        phi_list.append(phi_i)

    return phi_list


def pengrob_ki_list(
    pabs: float, tabs: float, ci_list: list, xi_list: list, yi_list: list, prop_dict: dict, bini_dict: dict
) -> list:
    """Peng Robinson Equilibrium Constants

    Calculate the peng robinson equilibrium constants for a mixture.

    Args:
        pabs (float): Absolute Evaluation Pressure, psia
        tabs (float): Absolute Evaluation Temp, rankine
        ci_list (list): List of String Components
        xi_list (list): Liquid Phase Molar Fractions
        yi_list (list): Vapor Phase Molar Fractions
        prop_dict (dict): Properties Dictionary
        bini_dict (dict): Binary Interaction Parameter Dictionary

    Returns:
        ki_list (float): Peng Robinson Equilibrium Constants
    """
    phi_vap_list = pengrob_fugco_list(pabs, tabs, ci_list, yi_list, prop_dict, bini_dict, True)
    phi_liq_list = pengrob_fugco_list(pabs, tabs, ci_list, xi_list, prop_dict, bini_dict, False)

    # calculate out a new ki, the book says to use log rules? can you get zero phi's?
    ki_list = [phi_liq / phi_vap for phi_vap, phi_liq in zip(phi_vap_list, phi_liq_list)]
    return ki_list
