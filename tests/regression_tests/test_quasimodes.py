from .regression import RegressionTest
import pytest
import numpy as np


class TestQuasimodesQR(RegressionTest):
    name = "quasimodes k2=1 k3=0.05"
    filename = "quasimodes_QR_k2_1_k3_0.05"
    equilibrium = "resonant_absorption"
    geometry = "Cartesian"

    parameters = {
        "k2": 1.0,
        "k3": 0.05,
        "p1": 0.9,
        "p2": 0.1,
        "r0": 0.1,
        "cte_T0": 0.0,
        "cte_B02": 0.0,
        "cte_B03": 1.0,
    }
    physics_settings = {
        "resistivity": True,
        "use_fixed_resistivity": True,
        "fixed_eta_value": 1e-4,
    }

    spectrum_limits = [
        {"xlim": (-60, 60), "ylim": (-8, 0.5)},
        {"xlim": (-0.5, 0.5), "ylim": (-1, 0.5)},
        {"xlim": (-0.18, 0.18), "ylim": (-0.15, 0.03)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    def test_eta_value(self, ds_test):
        assert np.all(
            ds_test.equilibria.get("eta")
            == pytest.approx(self.physics_settings["fixed_eta_value"])
        )
