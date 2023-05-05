import numpy as np
import pytest

from .regression import RegressionTest


class TestHarrisSheetQR(RegressionTest):
    name = "hall harris sheet k2=0.155 k3=0.01"
    filename = "hall_harris_sheet_QR_k2_0.155_k3_0.01"
    equilibrium = "harris_sheet"
    geometry = "Cartesian"
    x_start = -15
    x_end = 15

    parameters = {
        "k2": 0.155,
        "k3": 0.01,
        "alpha": 1,
        "cte_rho0": 1.0,
        "cte_B02": 1.0,
        "cte_B03": 5.0,
        "eq_bool": False,
    }
    physics_settings = {
        "resistivity": True,
        "fixed_resistivity_value": 1e-4,
        "hall_mhd": True,
        "electron_fraction": 0.5,
        "unit_density": 1.7e-14,
        "unit_magneticfield": 10,
        "unit_length": 7.534209349981049e-9,
        "incompressible": True,
    }

    spectrum_limits = [
        {"xlim": (-3.5e4, 3.5e4), "ylim": (-0.03, 0.013)},
        {"xlim": (-750, 750), "ylim": (-0.03, 0.013)},
        {"xlim": (-10, 10), "ylim": (-0.1, 0.05)},
        {"xlim": (-2.5, 2.5), "ylim": (-0.03, 0.03)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    def test_resistivity(self, ds_test):
        assert np.all(
            ds_test.equilibria.get("eta")
            == pytest.approx(self.physics_settings["fixed_resistivity_value"])
        )

    def test_units(self, ds_test):
        assert ds_test.cgs
        for val in ("density", "magneticfield", "length"):
            assert ds_test.units.get(f"unit_{val}") == pytest.approx(
                self.physics_settings[f"unit_{val}"]
            )

    def test_incompressible(self, ds_test, ds_base):
        assert ds_test.gamma == pytest.approx(ds_base.gamma)
        assert ds_test.gamma > 1e10

    def test_electron_fraction(self, ds_test, ds_base):
        key = "electronfraction"
        assert (
            ds_test.parameters[key]
            == pytest.approx(ds_base.parameters[key])
            == pytest.approx(0.5)
        )

    def test_electron_inertia(self, ds_test, ds_base):
        assert np.allclose(ds_test.equilibria["inertia"], 0)
        assert np.allclose(ds_base.equilibria["inertia"], 0)
