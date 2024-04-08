"""
What needs to be done?
(1)	Choose any overall composition for your selected mixture.
(2)	Calculate two saturation pressures (bubble and dew point) at two different temperatures using the Peng-Robinson EOS.
(3)	Finally, compare your calculated saturation pressures on the phase envelope for your respective mixtures.
    The phase envelope can be easily generated using the MI PVT software and the text output can be imported to Excel.

Useful hints/tips
(1)	For simplicity all binary interaction parameters can be set to zero.
(2)	A-priori generation of phase envelope using MI PVT will give you a clue
    if the selected temperatures will be bubble or dew points.
(3)	“Example calculation file in Excel” that is in the same folder as the final project has examples
    of saturation pressure and flash calculations that can be used/modified to accomplish your calculations.
"""

from num_methods import mix_comp_table
from overall import bubblepoint_pressure, comp_verify, dewpoint_pressure, phase_comp
from proptables.bini_vals import bini_dict
from proptables.crit_vals import prop_dict

comp = "c1"
# print(prop_dict[comp])

# print(f"Binary Interaction Parameter for C1 and nC5 is {bini_dict['c1']['nc5']}")

prac_comp = {"c3": 0.6, "nc4": 0.3, "nc5": 0.1}

# add some code to double check that it all sums to one, and it has recognized components

py_bub = bubblepoint_pressure(150, prac_comp, prop_dict, bini_dict)
print(f"\nPython BubblePoint: {round(py_bub,2)} psig\n")

py_dew = dewpoint_pressure(150, prac_comp, prop_dict, bini_dict)
print(f"Python DewPoint: {round(py_dew,2)} psig\n")

xi_list, yi_list = phase_comp(175, 150, prac_comp, prop_dict, bini_dict)
print(mix_comp_table(prac_comp, xi_list, yi_list))

test_comp = {"h2s": 0.1, "c3": 0.4, "nc4": 0.3, "nc5": 0.1, "nc8": 0.1}
comp_verify(test_comp, prop_dict, bini_dict)

xi_list, yi_list = phase_comp(175, 150, test_comp, prop_dict, bini_dict)
print(mix_comp_table(test_comp, xi_list, yi_list))
