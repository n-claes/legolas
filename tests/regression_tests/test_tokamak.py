tokamak_setup = {
    "name": "tokamak",
    "config": {
        "geometry": "cylindrical",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {
            "k2": -2.0,
            "k3": 0.2,
            "j0": (2.0 * 0.2) / 1.95,
            "cte_rho0": 1.0,
            "cte_B03": 1.0,
        },
        "equilibrium_type": "constant_current_tokamak",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "image_limits": [
        {"xlims": (-350, 350), "ylims": (-0.01, 0.01), "RMS_TOLERANCE": 2.2},
        {"xlims": (-0.01, 0.01), "ylims": (-0.01, 0.01), "RMS_TOLERANCE": 3.81},
    ],
}
