rti_thetapinch_hd_setup = {
    "name": "RTI_theta_pinch_HD",
    "config": {
        "gridpoints": 51,
        "parameters": {
            "k2": 1.0,
            "k3": 0.0,
            "cte_rho0": 1.0,
            "alpha": 2.0,
            "delta": 1 / 6,
            "r0": 0.0,
        },
        "flow": True,
        "equilibrium_type": "RTI_theta_pinch",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    }
}
