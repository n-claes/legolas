from .regression import RegressionTest
import pytest
import numpy as np


class TestRTI_KHI_QR(RegressionTest):
    name = "rayleigh-taylor kelvin-helmholtz k2=0 k3=1"
    filename = "rti_khi_QR_k2_0_k3_1"
    equilibrium = "RTI_KHI"
    geometry = "Cartesian"

    parameters = {
        "k2": 0.0,
        "k3": 1.0,
        "cte_rho0": 1.0,
        "cte_p0": 1000.0,
        "delta": -5.0,
        "g": 100.0,
        "alpha": -np.pi,
        "theta": 0.0,
        "p1": 1.0,
        "p2": 2.0,
        "p3": 1.0,
        "p4": 0.5 * np.pi,
        "tau": 4.0,
    }
    physics_settings = {"flow": True, "external_gravity": True}

    spectrum_limits = [
        {"xlim": (-1.3e4, 1.3e4), "ylim": (-4, 4)},
        {"xlim": (-350, 350), "ylim": (-4, 4)},
        {"xlim": (-5, 5), "ylim": (-4, 4)},
        {"xlim": (-2, 3.5), "ylim": (-4, 4)},
        {"xlim": (-1, 3), "ylim": (-0.6, 0.6)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    def test_gravity_value(self, ds_test):
        assert np.all(
            ds_test.equilibria.get("grav") == pytest.approx(self.parameters["g"])
        )
