from .regression import RegressionTest
import pytest


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
        "use_fixed_resistivity": True,
        "fixed_eta_value": 1e-4,
        "hall_mhd": True,
        "hall_substitution": True,
        "electron_fraction": 0.5,
        "cgs_units": True,
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
