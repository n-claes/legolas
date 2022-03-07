import numpy as np

NB_RUNS = 18
_kvals = np.linspace(0, np.sqrt(250), NB_RUNS)

gravito_hd_setup = {
    "name": "gravito_HD_multi",
    "config": {
        "equilibrium_type": "gravito_acoustic",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": _kvals,
            "k3": _kvals,
            "cte_p0": 1,
            "g": 0.5,
            "alpha": 20.42,
        },
        "external_gravity": True,
        "basename_datfile": "gravito_acoustic",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
    "xdata": 2 * _kvals**2,  # k0**2
    "y_scaling": 3 * 40.84 / 5,  # 1 / cs**2
    "use_squared_omega": True,
    "limits": {"xlims": (0, 550), "ylims": (0, 550)},
}
