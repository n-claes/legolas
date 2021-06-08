import pytest
import numpy as np

rti_setup = {
    "name": "rayleigh_taylor",
    "config": {
        "gridpoints": 51,
        "parameters": {
            "k2": 0.0,
            "k3": 1.0,
            "cte_rho0": 1.0,
            "cte_p0": 1000.0,
            "delta": -5.0,
            "g": 15.0,
            "alpha": 0.0,
            "theta": 0.35 * np.pi,
            "p1": 0.2,
            "p2": 0.6,
            "p3": 0.0,
            "p4": -0.35 * np.pi,
            "tau": 0.0,
        },
        "flow": True,
        "external_gravity": True,
        "equilibrium_type": "rayleigh_taylor",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "image_limits": [
        {"xlims": (-1.3e4, 1.3e4), "ylims": (-1.6, 1.6)},
        {"xlims": (-500, 500), "ylims": (-1.6, 1.6)},
        {"xlims": (-10, 15), "ylims": (-1.6, 1.6)},
        {"xlims": (0.01, 0.23), "ylims": (-1.6, 1.6)},
    ],
}
parametrisation = dict(
    argnames="setup",
    argvalues=[rti_setup],
    ids=[rti_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_gravity_value(ds_test, setup):
    g_value = 15.0
    assert setup["config"]["parameters"]["g"] == pytest.approx(g_value)
    assert np.all(ds_test.equilibria.get("grav") == pytest.approx(g_value))
