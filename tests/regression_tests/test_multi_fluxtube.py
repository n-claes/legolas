from .regression import MultiRegressionTest
import numpy as np


class TestFluxtubePhotosphericMulti(MultiRegressionTest):
    name = "fluxtube photospheric multi"
    filename = "multi_fluxtube_photospheric"
    equilibrium = "photospheric_flux_tube"
    geometry = "cylindrical"
    number_of_runs = 18
    gridpoints = 31

    x_start = 0
    x_end = 10
    r0 = 1.0

    k3vals = np.linspace(0.1, 6.2, number_of_runs)
    cs_min = 1.29099445

    parameters = {"k2": 0.0, "k3": k3vals, "cte_rho0": 1.0, "cte_p0": 1.0, "r0": r0}

    multispectrum_settings = {
        "xdata": k3vals * r0,
        "y_scaling": 1 / (k3vals * r0 * cs_min),
        "use_squared_omega": False,
        "xlim": (0, 6.3),
        "ylim": (0.84, 1.55),
    }


class TestFluxtubeCoronalMulti(MultiRegressionTest):
    name = "fluxtube coronal multi"
    filename = "multi_fluxtube_coronal"
    equilibrium = "coronal_flux_tube"
    geometry = "cylindrical"
    number_of_runs = 18
    gridpoints = 31

    x_start = 0
    x_end = 10
    r0 = 1.0

    k3vals = np.linspace(0.1, 6.2, number_of_runs)
    cs_max = 1.29099445

    parameters = {"k2": 0.0, "k3": k3vals, "cte_rho0": 1.0, "cte_p0": 1.0, "r0": r0}

    multispectrum_settings = {
        "xdata": k3vals * r0,
        "y_scaling": 1 / (k3vals * r0 * cs_max),
        "use_squared_omega": False,
        "xlim": (0, 6.3),
        "ylim": (0.84, 5.1),
    }
