import pytest
import numpy as np

gold_hoyle_setup = {
    "name": "gold_hoyle",
    "config": {
        "geometry": "cylindrical",
        "x_start": 0,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {
            "k2": 1.0,
            "k3": 1.0,
            "cte_rho0": 1.0,
            "cte_T0": 0.001,
            "alpha": 20.0,
        },
        "equilibrium_type": "gold_hoyle",
        "radiative_cooling": True,
        "cooling_curve": "rosner",
        "thermal_conduction": True,
        "use_fixed_tc_perp": True,
        "fixed_tc_perp_value": 0,
        "cgs_units": True,
        "unit_density": 1.6727e-15,
        "unit_magneticfield": 22.5,
        "unit_length": 1.0e10,
        "mean_molecular_weight": 1.0,
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "image_limits": [
        {"xlims": (-275, 275), "ylims": (-0.05, 1.2)},
        {"xlims": (-0.1, 0.1), "ylims": (0.92, 1)},
        {"xlims": (-150, 150), "ylims": (-0.005, 0.02)},
        {"xlims": (-20, 20), "ylims": (-0.0025, 0.018)},
        {"xlims": (-1.5, 1.5), "ylims": (-0.001, 0.018)},
    ],
}
parametrisation = dict(
    argnames="setup",
    argvalues=[gold_hoyle_setup],
    ids=[gold_hoyle_setup["name"]],
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


@pytest.mark.parametrize(**parametrisation)
def test_temperature(ds_test, setup):
    # relative tolerance of 0.2%
    assert np.all(
        ds_test.equilibria.get("T0") * ds_test.units.get("unit_temperature")
        == pytest.approx(2.9e5, rel=0.0065)
    )
