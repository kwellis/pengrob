import numpy as np

from num_methods import mix_comp_table
from overall import bubblepoint_pressure, comp_verify, dewpoint_pressure, phase_comp
from proptables.bini_vals import bini_dict
from proptables.crit_vals import prop_dict

# bad_comp = {"c3": 0.7, "nc4": 0.3, "nc5": 0.1} # does not sum to one
# bad_comp = {"h3": 0.6, "h4": 0.3, "nc5": 0.1} # components are not present
# comp_verify(bad_comp, prop_dict, bini_dict)

prac_comp = {"c3": 0.6, "nc4": 0.3, "nc5": 0.1}

# add some code to double check that it all sums to one, and it has recognized components
comp_verify(prac_comp, prop_dict, bini_dict)

py_bub = bubblepoint_pressure(100, prac_comp, prop_dict, bini_dict)
hy_bub = 124.9 - 14.7  # psia, hysys calculation, for composition
mi_bub = 134.28  # zero for binary interaction parameters
print(f"\nBubblePoints - MI PVT: {mi_bub} psi, Hysys: {hy_bub} psig, Python: {round(py_bub,2)} psig")

py_dew = dewpoint_pressure(100, prac_comp, prop_dict, bini_dict)
hy_dew = 67.3 - 14.7  # psia, hysys calculation for composition
mi_dew = np.nan
print(f"\nDewPoint - MI PVT: {mi_dew} psia, Hysys: {round(hy_dew,2)} psig, Python: {round(py_dew,2)} psig\n")

# need to build out a python check / verification for this
xi_list, yi_list = phase_comp(175, 150, prac_comp, prop_dict, bini_dict)
print(mix_comp_table(prac_comp, xi_list, yi_list))
