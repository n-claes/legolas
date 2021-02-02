import pytest
import numpy as np

interchange_modes_setup = {
    "name": "interchange_modes",
    "config": {
        "geometry": "Cartesian",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {
            "k2": np.pi,
            "k3": np.pi,
            "g": 0.5,
            "cte_p0": 0.25,
            "lambda": 0,
            "alpha": 20.0,
        },
        "equilibrium_type": "interchange_modes",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    }
}
parametrisation = dict(
    argnames="setup",
    argvalues=[interchange_modes_setup],
    ids=[interchange_modes_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_gravity_value(ds_test, setup):
    g_value = 0.5
    assert setup["config"]["parameters"]["g"] == pytest.approx(g_value)
    assert np.all(ds_test.equilibria.get("grav") == pytest.approx(g_value))
