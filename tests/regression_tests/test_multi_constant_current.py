import numpy as np

NB_RUNS = 24
j0 = (2.0 * 0.2) / np.linspace(1.9, 2.1, NB_RUNS)

constant_current_setup = {
    "name": "constant_current_multi",
    "config": {
        "equilibrium_type": "constant_current_tokamak",
        "gridpoints": 51,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": -2.0,
            "k3": 0.2,
            "j0": j0,
            "cte_rho0": 1.0,
            "cte_B03": 1.0,
        },
        "basename_datfile": "constant_current",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
    "xdata": 2 * 0.2 / j0,
    "symlog": 1e-8,
    "use_squared_omega": True,
    "limits": {"xlims": (1.88, 2.12), "ylims": (-1e-3, 1e6)},
}
