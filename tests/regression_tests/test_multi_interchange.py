import numpy as np

NB_RUNS = 24
ca_avg = 0.18257419

interchange_setup = {
    "name": "interchange_multi",
    "config": {
        "equilibrium_type": "interchange_modes",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": np.pi * np.sin(np.linspace(0, np.pi, NB_RUNS)),
            "k3": np.pi * np.cos(np.linspace(0, np.pi, NB_RUNS)),
            "g": 0.5,
            "cte_p0": 0.25,
            "lambda": 0.3,
            "alpha": 20.0,
        },
        "external_gravity": True,
        "basename_datfile": "interchange",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
    "xdata": np.linspace(0, np.pi, NB_RUNS) / np.pi,
    "y_scaling": 1 / ca_avg**2,
    "use_squared_omega": True,
    "limits": {"xlims": (-0.01, 1.01), "ylims": (-4.1, 14.4)},
}
