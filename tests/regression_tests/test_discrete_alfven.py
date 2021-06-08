import pytest
import numpy as np

discrete_alfven_setup = {
    "name": "discrete_alfven",
    "config": {
        "geometry": "cylindrical",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {"k2": 1.0, "k3": 0.05, "j0": 0.125, "delta": 0.2},
        "equilibrium_type": "discrete_alfven",
        "radiative_cooling": True,
        "cooling_curve": "rosner",
        "thermal_conduction": True,
        "use_fixed_tc_perp": True,
        "fixed_tc_perp_value": 0,
        "cgs_units": True,
        "unit_density": 1.5e-15,
        "unit_magneticfield": 50.0,
        "unit_length": 1.0e10,
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "image_limits": [
        {"xlims": (-700, 700), "ylims": (-0.9, 0.4)},
        {"xlims": (-100, 100), "ylims": (-0.9, 0.4)},
        {"xlims": (-0.2, 0.2), "ylims": (-0.9, 0.4)},
        {"xlims": (-0.01, 0.01), "ylims": (-0.004, 0.009), "RMS_TOLERANCE": 2.1},
    ],
}
parametrisation = dict(
    argnames="setup",
    argvalues=[discrete_alfven_setup],
    ids=[discrete_alfven_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
def test_conduction(ds_test, setup):
    assert np.all(
        ds_test.equilibria.get("kappa_perp") == setup["config"]["fixed_tc_perp_value"]
    )


@pytest.mark.parametrize(**parametrisation)
def test_units(ds_test, setup):
    assert ds_test.cgs
    for unit in ("unit_density", "unit_magneticfield", "unit_length"):
        assert ds_test.units.get(unit) == setup["config"][unit]
