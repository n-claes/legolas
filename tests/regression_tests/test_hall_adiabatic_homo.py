import numpy as np

hall_adiabatic_homo = {
    "name": "hall_adiabatic_homo",
    "config": {
        "geometry": "Cartesian",
        "x_start": 0,
        "x_end": 1000,
        "gridpoints": 51,
        "parameters": {
            "k2": round(np.pi * np.sin(np.pi / 6), 14),
            "k3": round(np.pi * np.cos(np.pi / 6), 14),
            "cte_rho0": 1.0,
            "cte_T0": 1.0,
            "cte_B02": 0.0,
            "cte_B03": 1.0,
        },
        "equilibrium_type": "adiabatic_homo",
        "hall_mhd": True,
        "electron_fraction": 0.5,
        "cgs_units": True,
        "unit_density": 1.7e-14,
        "unit_magneticfield": 10,
        "unit_length": 7.534209349981049e-9,
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": True,
    },
    "all_eigenvalues_real": True,
    "image_limits": [
        {"xlims": (-600, 600), "ylims": (-0.05, 0.05), "RMS_TOLERANCE": 2.5},
        {"xlims": (-50, 50), "ylims": (-0.05, 0.05)},
    ],
    "eigenfunctions": [
        {"eigenvalue": 0.7877457},
        {"eigenvalue": 4.0167199},
        {"eigenvalue": 9.4880654},
    ],
}
