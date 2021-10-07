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
    "all_eigenvalues_real": True,
    "image_limits": [
        {"xlims": (-600, 600), "ylims": (-0.05, 0.05)},
        {"xlims": (-50, 50), "ylims": (-0.05, 0.05)},
        {"xlims": (-0.5, 5), "ylims": (-0.05, 0.05)},
    ],
    "eigenfunctions": [
        {"eigenvalue": 2.67131},
        {"eigenvalue": 2.54724},
        {"eigenvalue": 2.51402},
        {"eigenvalue": 2.50119},
        {"eigenvalue": 2.49502},
    ],
}
