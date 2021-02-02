import pytest
import numpy as np

magnetothermal_setup = {
    "name": "magnetothermal",
    "config": {
        "geometry": "cylindrical",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {
            "k2": 0.0,
            "k3": 1.0,
            "cte_T0": 1.0,
        },
        "equilibrium_type": "magnetothermal_instabilities",
        "radiative_cooling": True,
        "cooling_curve": "rosner",
        "thermal_conduction": True,
        "use_fixed_tc_perp": True,
        "fixed_tc_perp_value": 0,
        "cgs_units": True,
        "unit_temperature": 2.6e6,
        "unit_magneticfield": 10.0,
        "unit_length": 1.0e8,
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    }
}
parametrisation = dict(
    argnames="setup",
    argvalues=[magnetothermal_setup],
    ids=[magnetothermal_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_conduction(ds_test, setup):
    assert np.all(
        ds_test.equilibria.get("kappa_perp") == setup["config"]["fixed_tc_perp_value"]
    )


@pytest.mark.parametrize(**parametrisation)
def test_units(ds_test, setup):
    assert ds_test.cgs
    for unit in ("unit_temperature", "unit_magneticfield", "unit_length"):
        assert ds_test.units.get(unit) == setup["config"][unit]
