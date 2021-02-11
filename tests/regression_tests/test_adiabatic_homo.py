import numpy as np

adiabatic_homo_setup = {
    "name": "adiabatic_homo",
    "config": {
        "geometry": "Cartesian",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {
            "k2": 0,
            "k3": np.pi,
            "cte_rho0": 1.0,
            "cte_T0": 1.0,
            "cte_B02": 0.0,
            "cte_B03": 1.0,
        },
        "equilibrium_type": "adiabatic_homo",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": True,
        "write_matrices": False,
    },
    "ev_guesses": [2.67131, 2.54724, 2.51402, 2.50119, 2.49502],
    "all_eigenvalues_real": True,
}
