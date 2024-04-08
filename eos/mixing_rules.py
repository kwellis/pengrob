def bini_aij(ai: float, aj: float, kij: float) -> float:
    """Binary Interaction of ai and aj

    Args:
        ai (float): PR or SRK a value for component i
        aj (float): PR or SRK a value for component j
        kij (float): Binary Interaction Parameter between Component i and j

    Returns:
        aij (float): Mixture PR or SRK a value
    """
    aij = (1 - kij) * (ai * aj) ** (1 / 2)
    return aij


def mix_a(ci_list: list, zi_list: list, ai_list: list, bini_dict: dict) -> float:
    """Mixture little a value

    Mixture value for little a, used for either peng robinson or srk

    Args:
        ci_list (list): Components in the Mixture, strings
        zi_list (list): Molar Fraction of the Components, floats
        ai_list (list): PR or SRK a parameters for each component
        bini_dict (dict): Dictionary of Binary Interaction Parameters

    Returns:
        mixa (float): Mixture Peng Robinson a Value
    """
    cj_list = ci_list.copy()
    zj_list = zi_list.copy()
    aj_list = ai_list.copy()

    at_list = []  # store all the values to be summed up

    # double for loop for double summation
    for ci, zi, ai in zip(ci_list, zi_list, ai_list):
        for cj, zj, aj in zip(cj_list, zj_list, aj_list):

            kij = bini_dict[ci][cj]  # binary interation parameter
            aij = bini_aij(ai, aj, kij)
            at = zi * zj * aij
            at_list.append(at)

    mixa = sum(at_list)
    return mixa


def mix_b(zi_list: list, bi_list: list) -> float:
    """Peng Robinson factor b for the mixture

    You could combine the functions for mixture functions of a and b if you wanted.
    Then it would output a float that represents the mixture for a and a mixture for b.

    Args:
        zi_list (list): Molar Fraction of the Components, floats
        bi_list (list): PR or SRK b parameters for each component

    Returns:
        mixb (float): Mixture Peng Robinson b Value
    """
    # update this code to mirror methods in other code?
    bt_list = []
    for zi, bi in zip(zi_list, bi_list):
        bt = zi * bi
        bt_list.append(bt)

    mixb = sum(bt_list)
    return mixb
