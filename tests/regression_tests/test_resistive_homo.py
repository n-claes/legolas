import pytest
import numpy as np

resistive_homo_setup = {
    "name": "resistive_homo",
    "config": {
        "geometry": "Cartesian",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 51,
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
    },
    "ev_guesses": [
        -0.4180167 - 0.0008883j,
        -0.4159505 - 0.0034422j,
        -0.4154890 - 0.0076967j,
        -0.4151935 - 0.0136528j,
        -0.4148061 - 0.0213107,
    ],
}
parametrisation = dict(
    argnames="setup",
    argvalues=[resistive_homo_setup],
    ids=[resistive_homo_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_eta_value(ds_test, setup):
    eta_value = 1e-3
    assert setup["config"]["fixed_eta_value"] == pytest.approx(eta_value)
    assert np.all(ds_test.equilibria.get("eta") == pytest.approx(eta_value))
