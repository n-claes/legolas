import numpy as np

NB_RUNS = 18
_kvals = np.linspace(0, np.sqrt(250), NB_RUNS)

gravito_mhd_setup = {
    "name": "gravito_MHD_multi",
    "config": {
        "equilibrium_type": "gravito_mhd",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": _kvals,
            "k3": _kvals,
            "cte_p0": 0.5,
            "g": 20,
            "alpha": 20,
        },
        "external_gravity": True,
        "basename_datfile": "gravito_mhd",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
    "xdata": 2 * _kvals**2,
    "y_scaling": 1,
    "use_squared_omega": True,
    "limits": {"xlims": (0, 550), "ylims": (0, 550)},
}
