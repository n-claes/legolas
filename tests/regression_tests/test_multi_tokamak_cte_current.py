from .regression import MultiRegressionTest
import numpy as np


class TestTokamakConstantCurrentMulti(MultiRegressionTest):
    name = "tokamak constant current multi"
    filename = "multi_tokamak_cte_current"
    equilibrium = "constant_current_tokamak"
    geometry = "cylindrical"
    number_of_runs = 24
    gridpoints = 51

    j0 = (2.0 * 0.2) / np.linspace(1.9, 2.1, number_of_runs)

    parameters = {"k2": -2.0, "k3": 0.2, "j0": j0, "cte_rho0": 1.0, "cte_B03": 1.0}

    multispectrum_settings = {
        "xdata": 2 * 0.2 / j0,
        "xlim": (1.88, 2.12),
        "ylim": (-1e-3, 1e6),
        "symlog": 1e-8,
    }
