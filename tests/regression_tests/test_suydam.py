suydam_setup = {
    "name": "suydam",
    "config": {
        "geometry": "cylindrical",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {
            "k2": 1.0,
            "k3": -1.2,
            "cte_rho0": 1,
            "cte_v02": 0,
            "cte_v03": 0.14,
            "cte_p0": 0.05,
            "p1": 0.1,
            "alpha": 2.0,
        },
        "flow": True,
        "equilibrium_type": "suydam_cluster",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "image_limits": [
        {"xlims": (-350, 350), "ylims": (-0.01, 0.01), "RMS_TOLERANCE": 2.25},
        {"xlims": (-1, 1), "ylims": (-0.01, 0.01)},
        {"xlims": (-0.2, -0.02), "ylims": (-0.01, 0.01), "RMS_TOLERANCE": 2.25},
        {"xlims": (-0.16, -0.08), "ylims": (-0.01, 0.01)},
    ]
}
