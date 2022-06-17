import numpy as np

NB_RUNS = 24
_kvals = np.logspace(-1, 4, NB_RUNS)

elecinertia_setup = {
    "name": "elecinertia_multi",
    "config": {
        "equilibrium_type": "adiabatic_homo",
        "gridpoints": 51,
        "number_of_runs": NB_RUNS,
        "geometry": "Cartesian",
        "x_start": 0,
        "x_end": 1000,
        "parameters": {
            "k2": _kvals * np.sin(np.pi/6),
            "k3": _kvals * np.cos(np.pi/6),
            "cte_rho0": 1,
            "cte_T0": 1,
            "cte_B02": 0,
            "cte_B03": 1
        },
        "hall_mhd": True,
        "hall_substitution": True,
        "electron_fraction": 0.5,
        "elec_inertia": True,
        "cgs_units": True,
        "unit_density": 1.7e-14,
        "unit_magneticfield": 10,
        "unit_length": 7.534209349981049e-9,
        "basename_datfile": "elecinertia",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
    "xdata": _kvals,
    "use_squared_omega": False,
    "xlog": True,
    "ylog": True,
    "limits": {"xlims": (1e-1, 1e4), "ylims": (1e-2, 1e4)},
}
