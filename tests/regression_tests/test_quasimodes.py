import pytest
import numpy as np

quasimodes_setup = {
    "name": "quasimodes",
    "config": {
        "geometry": "Cartesian",
        "x_start": 0.0,
        "x_end": 1.0,
        "gridpoints": 51,
        "parameters": {
            "k2": 1.0,
            "k3": 0.05,
            "p1": 0.9,
            "p2": 0.1,
            "r0": 0.1,
            "cte_T0": 0.0,
            "cte_B02": 0.0,
            "cte_B03": 1.0,
        },
        "equilibrium_type": "resonant_absorption",
        "resistivity": True,
        "use_fixed_resistivity": True,
        "fixed_eta_value": 10 ** (-4.0),
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "image_limits": [
        {"xlims": (-60, 60), "ylims": (-8, 0.5)},
        {"xlims": (-0.5, 0.5), "ylims": (-1, 0.5)},
        {"xlims": (-0.18, 0.18), "ylims": (-0.15, 0.03)},
    ],
}
parametrisation = dict(
    argnames="setup",
    argvalues=[quasimodes_setup],
    ids=[quasimodes_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_eta_value(ds_test, setup):
    eta_value = 1e-4
    assert setup["config"]["fixed_eta_value"] == pytest.approx(eta_value)
    assert np.all(ds_test.equilibria.get("eta") == pytest.approx(eta_value))
