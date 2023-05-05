import numpy as np

from .regression import MultiRegressionTest


class TestElectronInertiaMulti(MultiRegressionTest):
    name = "electron inertia multi"
    filename = "multi_elecinertia"
    equilibrium = "adiabatic_homo"
    geometry = "Cartesian"
    number_of_runs = 24
    x_start = 0
    x_end = 1000

    kvals = np.logspace(-1, 4, number_of_runs)

    parameters = {
        "k2": kvals * np.sin(np.pi / 6),
        "k3": kvals * np.cos(np.pi / 6),
        "cte_rho0": 1,
        "cte_T0": 1,
        "cte_B02": 0,
        "cte_B03": 1,
    }
    physics_settings = {
        "hall_mhd": True,
        "electron_fraction": 0.5,
        "elec_inertia": True,
        "unit_density": 1.7e-14,
        "unit_magneticfield": 10,
        "unit_length": 7.534209349981049e-9,
    }

    multispectrum_settings = {
        "xdata": kvals,
        "use_squared_omega": False,
        "xlog": True,
        "ylog": True,
        "xlim": (1e-1, 1e4),
        "ylim": (1e-2, 1e4),
    }
