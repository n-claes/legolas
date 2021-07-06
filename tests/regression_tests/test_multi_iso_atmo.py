import numpy as np

NB_RUNS = 12
NB_SCALEHEIGHTS = 15
SCALE_HEIGHT = 1.0
BETA2 = 0.5

temp_unit = 1e6
bfield_unit = 10.0
length_unit = 5e9
pres_unit = bfield_unit ** 2 / (4 * np.pi)
rho_unit = pres_unit * 1.672621777e-24 / (1.3806488e-16 * temp_unit)
vel_unit = bfield_unit / np.sqrt(4 * np.pi * rho_unit)
time_unit = length_unit / vel_unit

g = 2.74e4 / (length_unit / time_unit ** 2)
cte_T0 = g * SCALE_HEIGHT
gamma = 5 / 3

cs_max = 1.66319022

iso_atmo_beta_half_setup = {
    "name": "nye_thomas_beta_half",
    "config": {
        "x_start": 0,
        "x_end": NB_SCALEHEIGHTS * SCALE_HEIGHT,
        "equilibrium_type": "isothermal_atmosphere",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": 0,
            "k3": np.linspace(0, 5, NB_RUNS),
            "cte_rho0": 1,
            "cte_T0": cte_T0,
            "cte_B02": 0,
            "cte_B03": np.sqrt(gamma * cte_T0) * np.sqrt(BETA2),
            "g": g,
        },
        "external_gravity": True,
        "basename_datfile": "iso_atmo",
        "write_eigenfunctions": False,
        "logging_level": 0,
        "show_results": False,
        "unit_temperature": temp_unit,
        "unit_magneticfield": bfield_unit,
        "unit_length": length_unit,
    },
    "xdata": "k3",
    "y_scaling": SCALE_HEIGHT / cs_max,
    "use_squared_omega": False,
    "limits": {"xlims": (-0.2, 5.2), "ylims": (0, 10)},
}
