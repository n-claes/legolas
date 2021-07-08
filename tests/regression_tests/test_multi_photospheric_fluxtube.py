import numpy as np

NB_RUNS = 18
_k3vals = np.linspace(0.1, 6.2, NB_RUNS)
_r0 = 1.0
cs_min = 1.29099445

photospheric_tube_setup = {
    "name": "fluxtube_photo_series",
    "config": {
        "equilibrium_type": "photospheric_flux_tube",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": 0.0,
            "k3": _k3vals,
            "cte_rho0": 1.0,
            "cte_p0": 1.0,
            "r0": _r0,
        },
        "basename_datfile": "photo_tube",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
    "xdata": _k3vals * _r0,
    "y_scaling": 1 / (_k3vals * _r0 * cs_min),
    "use_squared_omega": False,
    "limits": {"xlims": (0, 6.3), "ylims": (0.84, 1.55)},
}
