mri_setup = {
    "name": "MRI_accretion",
    "config": {
        "gridpoints": 51,
        "parameters": {"k2": 0.0, "k3": 70.0, "beta": 100.0, "tau": 1.0, "nu": 0.1},
        "flow": True,
        "external_gravity": True,
        "equilibrium_type": "MRI_accretion",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": True,
        "write_matrices": False,
    },
    "ev_guesses": [
        -0.002028 + 0.627722j,
        -0.00186038 + 0.5804956j,
        -0.0017397 + 0.5443967j,
    ],
    "image_limits": [
        {"xlims": (-45, 45), "ylims": (-0.7, 0.7)},
        {"xlims": (-10, 10), "ylims": (-0.7, 0.7)},
        {"xlims": (-1.5, 1.5), "ylims": (-0.7, 0.7)},
        {"xlims": (-0.05, 0.05), "ylims": (-0.7, 0.7)},
        {"xlims": (-0.01, 0.01), "ylims": (-0.01, 0.7)},
    ]
}
