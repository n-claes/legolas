from .regression import MultiRegressionTest
import numpy as np


class TestNyeThomasMulti(MultiRegressionTest):
    name = "Nye Thomas isothermal atmosphere multi beta=0.5"
    filename = "multi_nye_thomas_iso_atmo_beta_half"
    equilibrium = "isothermal_atmosphere"
    geometry = "Cartesian"
    number_of_runs = 12
    gridpoints = 31

    nb_scaleheights = 15
    scaleheight = 1.0

    x_start = 0
    x_end = nb_scaleheights

    beta2 = 0.5
    temp_unit = 1e6
    bfield_unit = 10.0
    length_unit = 5e9
    pres_unit = bfield_unit**2 / (4 * np.pi)
    rho_unit = pres_unit * 1.672621777e-24 / (1.3806488e-16 * temp_unit)
    vel_unit = bfield_unit / np.sqrt(4 * np.pi * rho_unit)
    time_unit = length_unit / vel_unit

    g = 2.74e4 / (length_unit / time_unit**2)
    cte_T0 = g * scaleheight
    gamma = 5 / 3
    cs_max = 1.66319022

    parameters = {
        "k2": 0,
        "k3": np.linspace(0, 5, number_of_runs),
        "cte_rho0": 1,
        "cte_T0": cte_T0,
        "cte_B02": 0,
        "cte_B03": np.sqrt(gamma * cte_T0) * np.sqrt(beta2),
        "g": g,
    }
    physics_settings = {
        "external_gravity": True,
        "unit_temperature": temp_unit,
        "unit_magneticfield": bfield_unit,
        "unit_length": length_unit,
    }

    multispectrum_settings = {
        "xdata": "k3",
        "y_scaling": scaleheight / cs_max,
        "use_squared_omega": False,
        "xlim": (-0.2, 5.2),
        "ylim": (0, 10),
    }
