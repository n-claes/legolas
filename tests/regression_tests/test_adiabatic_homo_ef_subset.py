import numpy as np

adiabatic_ef_subset_setup = {
    "name": "adiabatic_homo_ef_subset",
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
        "write_eigenfunction_subset": True,
        "eigenfunction_subset_center": 20 + 1j,
        "eigenfunction_subset_radius": 15,
    },
    "all_eigenvalues_real": True,
    "eigenfunctions": [
        {"eigenvalue": 11.18509},
        {"eigenvalue": 16.02714},
        {"eigenvalue": 21.00395},
    ],
}
