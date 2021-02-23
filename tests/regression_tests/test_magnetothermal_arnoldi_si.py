import pytest
import numpy as np

magneto_arnoldi_si_setup = {
    "name": "magnetothermal_arnoldi_si",
    "config": {
        "geometry": "cylindrical",
        "x_start": 0.01,
        "x_end": 1,
        "gridpoints": 50,
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
        "write_eigenfunctions": True,
        "write_matrices": False,
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 100,
        "which_eigenvalues": "LM",
        "sigma": 0.1 + 0.05j,
    },
    "ev_guesses": [
        2.799e-2j,
        2.022e-2 + 4.418e-2j,
    ],
}
parametrisation = dict(
    argnames="setup",
    argvalues=[magneto_arnoldi_si_setup],
    ids=[magneto_arnoldi_si_setup["name"]],
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
