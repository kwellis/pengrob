import math

import numpy as np

import eos.eos_start as es
import eos.peng_robinson as pr
import num_methods as nm
import rachford_rice as rr


def comp_verify(comp_dict: dict, prop_dict: dict, bini_dict: dict) -> None:
    """Composition Verification

    Verify the molar composition sums to one and all the components are present

    Args:
        comp_dict (dict): Mixture Molar Composition
        prop_dict (dict): Property Table for Lookup
        bini_dict (dict): Binary Interaction Table

    Returns:
        None if no errors are present.
    """
    zi_list = comp_dict.values()
    zi_tot = sum(zi_list)
    zi_err = 1 - zi_tot
    err_tot = 1e-3  # how much off can the values be from adding to one

    if abs(zi_err) > err_tot:
        raise ValueError(f"Molar fractions do not sum to one, Error is {zi_err: .3E}")

    miss_prop = comp_dict.keys() - prop_dict.keys()
    miss_bini = comp_dict.keys() - bini_dict.keys()

    if bool(miss_prop) is True:
        raise KeyError(f"{miss_prop} are not defined in property table")

    if bool(miss_bini) is True:
        raise KeyError(f"{miss_bini} are not defined in binary interaction table")

    return None


def bubblepoint_pressure(teval: float, comp_dict: dict, prop_dict: dict, bini_dict: dict) -> float:
    """Peng Robinson Bubble Point Pressure

    Args:
        teval (float): Evaluation Temperature, deg F
        comp_dict (dict): Mixture Molar Composition
        prop_dict (dict): Property Table for Lookup
        bini_dict (dict): Binary Interaction Table

    Returns:
        pbub (float): Bubble Point Pressure, psig
    """
    comp_list = list(comp_dict.keys())
    zi_list = list(comp_dict.values())
    xi_list = zi_list.copy()

    tabs = teval + 459.67
    plist = [es.bubblepoint_guess(tabs, comp_dict, prop_dict)]  # starting / guess pressure

    ki_list = es.wilson_ki_list(plist[-1], tabs, comp_list, zi_list, prop_dict, bini_dict)
    yi_list = [rr.vapor_frac(zi, ki, 0) for zi, ki in zip(zi_list, ki_list)]

    ki_list = pr.pengrob_ki_list(plist[-1], tabs, comp_list, xi_list, yi_list, prop_dict, bini_dict)
    yi_list = [rr.vapor_frac(zi, ki, 0) for zi, ki in zip(zi_list, ki_list)]
    yi_tot_list = [sum(yi_list)]  # store this, do secant on this value

    # rough approximation for going up or down, before secant method takes over
    if yi_tot_list[0] > 0:
        plist.append(plist[0] + 50)
    else:
        plist.append(plist[0] - 50)

    pdiff = 0.001  # how much the iteration needs to change
    while abs(plist[-2] - plist[-1]) > pdiff:
        ki_list = pr.pengrob_ki_list(plist[-1], tabs, comp_list, xi_list, yi_list, prop_dict, bini_dict)
        yi_list = [rr.vapor_frac(zi, ki, 0) for zi, ki in zip(zi_list, ki_list)]
        yi_tot_list.append(sum(yi_list))

        plist.append(nm.psi_secant(plist[-2], plist[-1], yi_tot_list[-2], yi_tot_list[-1]))

    return plist[-1] - 14.7


def dewpoint_pressure(teval: float, comp_dict: dict, prop_dict: dict, bini_dict: dict) -> float:
    """Peng Robinson Dew Point Pressure

    Args:
        teval (float): Evaluation Temperature, deg F
        comp_dict (dict): Mixture Molar Composition
        prop_dict (dict): Property Table for Lookup
        bini_dict (dict): Binary Interaction Table

    Returns:
        pdew (float): Dew Point Pressure, psig
    """
    comp_list = list(comp_dict.keys())
    zi_list = list(comp_dict.values())
    yi_list = zi_list.copy()

    tabs = teval + 459.67
    plist = [es.dewpoint_guess(tabs, comp_dict, prop_dict)]  # starting / guess pressure

    ki_list = es.wilson_ki_list(plist[-1], tabs, comp_list, zi_list, prop_dict, bini_dict)
    xi_list = [rr.liquid_frac(zi, ki, 1) for zi, ki in zip(zi_list, ki_list)]

    ki_list = pr.pengrob_ki_list(plist[-1], tabs, comp_list, xi_list, yi_list, prop_dict, bini_dict)
    xi_list = [rr.liquid_frac(zi, ki, 1) for zi, ki in zip(zi_list, ki_list)]
    xi_tot_list = [sum(xi_list)]

    # rough approximation for going up or down, before secant method takes over
    if xi_tot_list[0] > 0:
        plist.append(plist[0] + 50)
    else:
        plist.append(plist[0] - 50)

    pdiff = 0.001  # how much the iteration needs to change
    while abs(plist[-2] - plist[-1]) > pdiff:
        ki_list = pr.pengrob_ki_list(plist[-1], tabs, comp_list, xi_list, yi_list, prop_dict, bini_dict)
        xi_list = [rr.liquid_frac(zi, ki, 1) for zi, ki in zip(zi_list, ki_list)]
        xi_tot_list.append(sum(xi_list))

        plist.append(nm.psi_secant(plist[-2], plist[-1], xi_tot_list[-2], xi_tot_list[-1]))

    return plist[-1] - 14.7


def phase_comp(peval: float, teval: float, comp_dict: dict, prop_dict: dict, bini_dict: dict) -> tuple[list, list]:
    """Peng Robinson Two Phase Composition

    Input a feed composition at a certain pressure and temperature.
    Output the composition of the xi, the liquid and yi, the vapor

    Args:
        peval (float): Evaluated Pressure, psig
        teval (float): Evaluated Temperature, deg F
        comp_dict (dict): Mixture Molar Composition
        prop_dict (dict): Property Table for Lookup
        bini_dict (dict): Binary Interaction Table

    Returns:
        xi_list (list): Liquid Molar Fraction Composition
        yi_list (list): Vapour Molar Fraction Composition
    """
    tabs = teval + 459.67
    pabs = peval + 14.7
    ci_list = list(comp_dict.keys())
    zi_list = list(comp_dict.values())

    beta = [0.5]  # vapor mole fraction starting point

    ki_list = es.wilson_ki_list(
        pabs,
        tabs,
        ci_list,
        zi_list,
        prop_dict,
        bini_dict,
    )

    # use wilson to calculate xi and yi fractions, move to peng rob ki
    xi_list = [rr.liquid_frac(zi, ki, beta[-1]) for zi, ki in zip(zi_list, ki_list)]
    yi_list = [rr.vapor_frac(zi, ki, beta[-1]) for zi, ki in zip(zi_list, ki_list)]
    ki_list = pr.pengrob_ki_list(pabs, tabs, ci_list, xi_list, yi_list, prop_dict, bini_dict)

    rrf_tot, rrd_tot = rr.rr_sum(zi_list, ki_list, beta[-1])
    beta.append(rr.rr_newton(rrf_tot, rrd_tot, beta[-1]))  # calculate next beta

    bdiff = 1e-5  # how much iteration needs to change
    while abs(beta[-2] - beta[-1]) > bdiff:

        xi_list = [rr.liquid_frac(zi, ki, beta[-1]) for zi, ki in zip(zi_list, ki_list)]
        yi_list = [rr.vapor_frac(zi, ki, beta[-1]) for zi, ki in zip(zi_list, ki_list)]
        ki_list = pr.pengrob_ki_list(pabs, tabs, ci_list, xi_list, yi_list, prop_dict, bini_dict)

        rrf_tot, rrd_tot = rr.rr_sum(zi_list, ki_list, beta[-1])
        beta.append(rr.rr_newton(rrf_tot, rrd_tot, beta[-1]))  # calculate next beta

    return xi_list, yi_list
