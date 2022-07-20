from .regression import MultiRegressionTest
import numpy as np


class TestInterchangeModesMulti(MultiRegressionTest):
    name = "interchange modes multi"
    filename = "multi_interchange_modes"
    equilibrium = "interchange_modes"
    geometry = "Cartesian"
    number_of_runs = 24
    gridpoints = 31

    k2vals = np.pi * np.sin(np.linspace(0, np.pi, number_of_runs))
    k3vals = np.pi * np.cos(np.linspace(0, np.pi, number_of_runs))

    parameters = {
        "k2": k2vals,
        "k3": k3vals,
        "g": 0.5,
        "cte_p0": 0.25,
        "lambda": 0.3,
        "alpha": 20.0,
    }
    physics_settings = {"external_gravity": True}

    multispectrum_settings = {
        "xdata": np.linspace(0, np.pi, number_of_runs) / np.pi,
        "y_scaling": 1 / 0.18257419**2,  # 1 / cA**2
        "xlim": (-0.01, 1.01),
        "ylim": (-4.1, 14.4),
    }
