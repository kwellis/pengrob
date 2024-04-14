"""Visual Graphs and Plots

Used for creating visuals to be seen and verify results are correct
"""

import matplotlib.pyplot as plt
import pandas as pd

tern = pd.read_excel("data/hysys_vals.xlsx", sheet_name="ternary")
lift = pd.read_excel("data/hysys_vals.xlsx", sheet_name="lift_gas")


def hysys_plot(
    py_pres: list[float],
    py_temp: list[float],
    py_desc: list[str],
    hy_pres: list[float],
    hy_temp: list[float],
    hy_desc: list[str],
    title: str,
) -> None:
    """Hysys Plot and Visualization

    Create a plot to compare the python data to the hysys data

    Args:
        py_press
        py_temp
        py_desc
        hy_pres
        hy_temp
        hy_desc
        title

    Returns:
        Graphs (None)
    """
    for pres, temp, desc in zip(py_pres, py_temp, py_desc):
        plt.plot(temp, pres, marker="o", label="Py-" + desc.capitalize())

    dew_idx = hy_desc == "dew"
    plt.plot(hy_temp[dew_idx], hy_pres[dew_idx], label="Hy-Dew")

    bub_idx = hy_desc == "bub"
    plt.plot(hy_temp[bub_idx], hy_pres[bub_idx], label="Hy-Bub")

    plt.title(title)
    plt.xlabel("Temperature, Deg F")
    plt.ylabel("Pressure, PSIG")
    plt.legend()
    plt.show()
