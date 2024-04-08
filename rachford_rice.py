"""Rachford and Rice Objectives Functions and Derivatives

Solve the vapor molar fraction for flash equations, referred to as beta in Pedersen Book.
https://www.e-education.psu.edu/png520/m13_p2.html
"""


def rr_func(zi: float, Ki: float, beta: float) -> float:
    """Component Rachford Rice Function

    Args:
        zi (float): Feed Mixture Molar Fraction, decimal
        Ki (float): Equilibrium Ratio, fug coeff liquid / fug coeff vapor
        beta (float): Vapor Mole Fraction, Total Mixture

    Return:
        rrfi (float): Rachford Rice Function for single component
    """
    rrfi = zi * (Ki - 1) / (1 + beta * (Ki - 1))
    return rrfi


def rr_derv(zi: float, Ki: float, beta: float) -> float:
    """Component Rachford Rice Derivative

    Args:
        zi (float): Feed Mixture Molar Fraction, decimal
        Ki (float): Equilibrium Ratio, fug coeff liquid / fug coeff vapor
        beta (float): Vapor Mole Fraction, Total Mixture

    Return:
        rrdi (float): Rachford Rice Derivative for single component
    """
    rrdi = -zi * (Ki - 1) ** 2 / (1 + beta * (Ki - 1)) ** 2
    return rrdi


def rr_sum(zi_list: list, Ki_list: list, beta: float) -> tuple[float, float]:
    """Rachford and Rice Equation and Derivative Summations

    Args:
        zi_list (list): Feed Mixture Molar Fractions
        Ki_list (list): Equilibrium Ratios of Components, fugco_liq / fugco_vap
        beta (float): Vapor Mole Fraction, Total Mixture

    Return:
        rrf_sum (float): Rachford Rice Function Summation
        rrd_sum (float): Rachford Rice Derivative Summation
    """
    rrf_list = []
    rrd_list = []
    for zi, Ki in zip(zi_list, Ki_list):
        rrf_list.append(rr_func(zi, Ki, beta))
        rrd_list.append(rr_derv(zi, Ki, beta))
    return sum(rrf_list), sum(rrd_list)


def rr_newton(rrf_sum: float, rrd_sum: float, beta: float) -> float:
    """Rachford and Rice Newton Beta Calculation

    Args:
        rrf_sum (float): Rachford Rice Function Summation
        rrd_sum (float): Rachford Rice Derivative Summation
        beta (float): Vapor Mole Fraction, Total Mixture

    Return:
        beta_nxt (float): Next Vapor Mole Fraction
    """
    beta_nxt = beta - rrf_sum / rrd_sum
    return beta_nxt


def rr_subcool(zi_list: list, Ki_list: list) -> bool:
    """Sub Cooled Liquid Check

    Equation for testing for the presence of a subcooled liquid. If
    the mixture is a subcooled liquid it means that beta (vapor) is zero
    and the mixture is totally liquid, with no gas. Pg 252 Michelsen 2008

    Args:
        zi_list (list): Feed Mixture Molar Fractions
        Ki_list (list): Equilibrium Ratios of Components, fugco_liq / fugco_vap

    Return:
        subcool (bool): True - Mixture is totally liquid, False - Mixture has vapor
    """
    rrf, rrd = rr_sum(zi_list, Ki_list, 0)
    if rrf > 0:
        subcool = False
    else:
        subcool = True
    return subcool


def rr_superheat(zi_list: list, Ki_list: list) -> bool:
    """Super Heated Vapor Check

    Equation for testing for the presence of a super heated vapor. If
    the mixture is a super heated vapor it means that beta (vapor) is one
    and the mixture is totally vapor, with no liquid. Pg 252 Michelsen 2008

    Args:
        zi_list (list): Feed Mixture Molar Fractions
        Ki_list (list): Equilibrium Ratios of Components, fugco_liq / fugco_vap

    Return:
        supheat (bool): True - Mixture is totally vapor, False - Mixture has liquid
    """
    rrf, rrd = rr_sum(zi_list, Ki_list, 1)
    if rrf < 0:
        supheat = False
    else:
        supheat = True
    return supheat


def vapor_frac(zi: float, Ki: float, beta: float) -> float:
    """Vapor Fraction yi of Component

    Calculate vapor fraction of single component.

    Args:
        zi (float): Feed Molar Fraction for single component, decimal
        Ki (float): Equilibrium Ratio, fug coeff liquid / fug coeff vapor
        beta (float): Vapor Mole Fraction for total mixture, decimal

    Return:
        yi (float): Vapor Molar Fraction for single component, decimal
    """
    yi = Ki * zi / (1 - beta + beta * Ki)
    return yi


def liquid_frac(zi: float, Ki: float, beta: float) -> float:
    """Liquid Fraction xi of Component

    Calculate liquid fraction of single component.

    Args:
        zi (float): Feed Molar Fraction for single component, decimal
        Ki (float): Equilibrium Ratio, fug coeff liquid / fug coeff vapor
        beta (float): Vapor Mole Fraction for total mixture, decimal

    Return:
        xi (float): Liquid Molar Fraction for single component, decimal
    """
    xi = zi / (1 - beta + beta * Ki)
    return xi
