import pytest
import numpy as np

rti_khi_setup = {
    "name": "RTI_KHI",
    "config": {
        "gridpoints": 51,
        "parameters": {
            "k2": 0.0,
            "k3": 1.0,
            "cte_rho0": 1.0,
            "cte_p0": 1000.0,
            "delta": -5.0,
            "g": 100.0,
            "alpha": -np.pi,
            "theta": 0.0,
            "p1": 1.0,
            "p2": 2.0,
            "p3": 1.0,
            "p4": 0.5 * np.pi,
            "tau": 4.0,
        },
        "flow": True,
        "external_gravity": True,
        "equilibrium_type": "RTI_KHI",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    }
}
parametrisation = dict(
    argnames="setup",
    argvalues=[rti_khi_setup],
    ids=[rti_khi_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_gravity_value(ds_test, setup):
    g_value = 100.0
    assert setup["config"]["parameters"]["g"] == pytest.approx(g_value)
    assert np.all(ds_test.equilibria.get("grav") == pytest.approx(g_value))

