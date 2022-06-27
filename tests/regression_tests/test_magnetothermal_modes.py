from .regression import RegressionTest
import pytest
import numpy as np


class TestMagnetoThermalModesQR(RegressionTest):
    name = "magnetothermal modes k2=0 k3=1 QR"
    filename = "magnetothermal_QR_k2_0_k3_1"
    equilibrium = "magnetothermal_instabilities"
    geometry = "cylindrical"

    parameters = {"k2": 0.0, "k3": 1.0, "cte_T0": 1.0}
    physics_settings = {
        "radiative_cooling": True,
        "cooling_curve": "rosner",
        "thermal_conduction": True,
        "use_fixed_tc_perp": True,
        "fixed_tc_perp_value": 0,
        "unit_temperature": 2.6e6,
        "unit_magneticfield": 10.0,
        "unit_length": 1.0e8,
        "mean_molecular_weight": 1.0,
    }

    spectrum_limits = [
        {"xlim": (-650, 650), "ylim": (-0.15, 0.12)},
        {"xlim": (-40, 40), "ylim": (-0.15, 0.12)},
        {"xlim": (-5, 5), "ylim": (-0.15, 0.12)},
        {"xlim": (-0.1, 0.1), "ylim": (-0.15, 0.12)},
        {"xlim": (-0.025, 0.025), "ylim": (-0.015, 0.12), "RMS_TOLERANCE": 2.2},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    def test_conduction(self, ds_test):
        assert np.all(
            ds_test.equilibria.get("kappa_perp")
            == pytest.approx(self.physics_settings["fixed_tc_perp_value"])
        )

    def test_units(self, ds_test):
        assert ds_test.cgs
        for val in ("temperature", "magneticfield", "length"):
            assert (
                ds_test.units.get(f"unit_{val}") == self.physics_settings[f"unit_{val}"]
            )
