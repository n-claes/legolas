import pytest
import numpy as np

resistive_tearing_setup = {
    "name": "resistive_tearing",
    "config": {
        "geometry": "Cartesian",
        "x_start": -0.5,
        "x_end": 0.5,
        "gridpoints": 51,
        "parameters": {
            "k2": 0.49,
            "k3": 0.0,
            "alpha": 4.73884,
            "beta": 0.15,
            "cte_rho0": 1.0,
        },
        "equilibrium_type": "resistive_tearing",
        "resistivity": True,
        "use_fixed_resistivity": True,
        "fixed_eta_value": 10 ** (-4.0),
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
}
parametrisation = dict(
    argnames="setup",
    argvalues=[resistive_tearing_setup],
    ids=[resistive_tearing_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_eta_value(ds_test, setup):
    eta_value = 1e-4
    assert setup["config"]["fixed_eta_value"] == pytest.approx(eta_value)
    assert np.all(ds_test.equilibria.get("eta") == pytest.approx(eta_value))
