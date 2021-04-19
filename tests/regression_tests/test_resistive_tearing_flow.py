import pytest
import numpy as np

resistive_tearing_flow_setup = {
    "name": "resistive_tearing_flow",
    "config": {
        "geometry": "Cartesian",
        "x_start": -0.5,
        "x_end": 0.5,
        "gridpoints": 51,
        "parameters": {
            "k2": 1.5,
            "k3": 0.0,
            "alpha": 4.73884,
            "beta": 0.15,
            "cte_rho0": 1.0,
        },
        "equilibrium_type": "resistive_tearing_flow",
        "flow": True,
        "resistivity": True,
        "use_fixed_resistivity": True,
        "fixed_eta_value": 10 ** (-4.0),
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "image_limits": [
        {"xlims": (-375, 375), "ylims": (-11, 0.5)},
        {"xlims": (-30, 30), "ylims": (-3, 0.3)},
        {"xlims": (-4, 4), "ylims": (-2.2, 0.2)},
        {"xlims": (-1.6, 1.6), "ylims": (-1.6, 0.1)},
        {"xlims": (-0.6, 0.6), "ylims": (-0.38, 0.05)},
    ],
}
parametrisation = dict(
    argnames="setup",
    argvalues=[resistive_tearing_flow_setup],
    ids=[resistive_tearing_flow_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_eta_value(ds_test, setup):
    eta_value = 1e-4
    assert setup["config"]["fixed_eta_value"] == pytest.approx(eta_value)
    assert np.all(ds_test.equilibria.get("eta") == pytest.approx(eta_value))
