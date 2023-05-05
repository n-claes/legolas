import numpy as np
import pytest

from .regression import RegressionTest


class TestResistiveTearingQR(RegressionTest):
    name = "resistive tearing k2=0.49 k3=0"
    filename = "resistive_tearing_QR_k2_0.49_k3_0"
    equilibrium = "resistive_tearing"
    geometry = "Cartesian"
    x_start = -0.5
    x_end = 0.5

    parameters = {
        "k2": 0.49,
        "k3": 0.0,
        "alpha": 4.73884,
        "beta": 0.15,
        "cte_rho0": 1.0,
    }
    physics_settings = {
        "resistivity": True,
        "fixed_resistivity_value": 1e-4,
    }

    spectrum_limits = [
        {"xlim": (-375, 375), "ylim": (-20, 5)},
        {"xlim": (-30, 30), "ylim": (-15, 3)},
        {"xlim": (-10, 10), "ylim": (-7.5, 1)},
        {"xlim": (-1.2, 1.2), "ylim": (-1.1, 0.05)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    def test_eta_value(self, ds_test):
        assert np.all(
            ds_test.equilibria.get("eta")
            == pytest.approx(self.physics_settings["fixed_resistivity_value"])
        )


class TestResistiveTearingFlowQR(TestResistiveTearingQR):
    name = "resistive tearing flow k2=1.5 k3=0"
    filename = "resistive_tearing_flow_QR_k2_1.5_k3_0"
    equilibrium = "resistive_tearing_flow"

    parameters = {
        "k2": 1.5,
        "k3": 0.0,
        "alpha": 4.73884,
        "beta": 0.15,
        "cte_rho0": 1.0,
    }
    physics_settings = {
        "resistivity": True,
        "fixed_resistivity_value": 1e-4,
        "flow": True,
    }
    spectrum_limits = [
        {"xlim": (-375, 375), "ylim": (-11, 5)},
        {"xlim": (-30, 30), "ylim": (-3, 0.3)},
        {"xlim": (-4, 4), "ylim": (-2.2, 0.2)},
        {"xlim": (-1.6, 1.6), "ylim": (-1.6, 0.1)},
        {"xlim": (-0.6, 0.6), "ylim": (-0.38, 0.05)},
    ]
