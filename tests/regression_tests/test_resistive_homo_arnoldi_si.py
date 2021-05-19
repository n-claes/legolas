import pytest
import numpy as np

resistive_homo_arnoldi_si_setup = {
    "name": "resistive_homo_arnoldi_si",
    "config": {
        "geometry": "Cartesian",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 50,
        "parameters": {
            "k2": 0.0,
            "k3": 1.0,
            "beta": 0.25,
            "cte_rho0": 1.0,
            "cte_B02": 0.0,
            "cte_B03": 1.0,
        },
        "equilibrium_type": "resistive_homo",
        "resistivity": True,
        "use_fixed_resistivity": True,
        "fixed_eta_value": 0.001,
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": True,
        "write_matrices": False,
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 100,
        "which_eigenvalues": "LM",
        "sigma": 1.0 - 1.5j,
    },
    "image_limits": [
        {"xlims": (-0.6, 1.6), "ylims": (-3, 0.1), "RMS_TOLERANCE": 2.3},
        {"xlims": (-0.07, 1.04), "ylims": (-1.1, 0.1), "RMS_TOLERANCE": 2.55},
    ],
    "eigenfunctions": [
        {"eigenvalue": 4.180e-01 - 8.884e-04j},
        {"eigenvalue": 4.160e-01 - 3.442e-03j},
        {"eigenvalue": 4.155e-01 - 7.697e-03j},
        {"eigenvalue": 4.152e-01 - 1.365e-02j},
    ],
}
parametrisation = dict(
    argnames="setup",
    argvalues=[resistive_homo_arnoldi_si_setup],
    ids=[resistive_homo_arnoldi_si_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_eta_value(ds_test, setup):
    eta_value = 1e-3
    assert setup["config"]["fixed_eta_value"] == pytest.approx(eta_value)
    assert np.all(ds_test.equilibria.get("eta") == pytest.approx(eta_value))
