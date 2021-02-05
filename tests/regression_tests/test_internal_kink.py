internal_kink_setup = {
    "name": "internal_kink",
    "config": {
        "geometry": "cylindrical",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {
            "k2": 1.0,
            "k3": 0.16 * 5,
            "cte_rho0": 1.0,
            "cte_v03": 1.0,
            "cte_p0": 3.0,
            "alpha": 5.0,
        },
        "equilibrium_type": "internal_kink",
        "flow": True,
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": True,
        "write_matrices": False,
    },
    "ev_guesses": [0.470629 + 0.00607j, 0.470629 + 0.00607j],
}
