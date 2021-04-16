fluxtube_coronal_setup = {
    "name": "fluxtube_coronal",
    "config": {
        "geometry": "cylindrical",
        "x_start": 0,
        "x_end": 10,
        "mesh_accumulation": True,
        "gridpoints": 51,
        "parameters": {
            "k2": 0,
            "k3": 4.0,
            "cte_rho0": 1.0,
            "cte_p0": 1.0,
            "r0": 1.0,
        },
        "equilibrium_type": "coronal_flux_tube",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "all_eigenvalues_real": True,
    "image_limits": [
        {"xlims": (-6300, 6300), "ylims": (-0.5, 0.5)},
        {"xlims": (-1000, 1000), "ylims": (-0.5, 0.5)},
        {"xlims": (-50, 50), "ylims": (-0.5, 0.5), "RMS_TOLERANCE": 2.2},
    ],
}
