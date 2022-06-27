from .regression import RegressionTest
import pytest
import numpy as np


class TestRTI_QR(RegressionTest):
    name = "RTI k2=1 k3=0.05"
    filename = "rti_QR_k2_1_k3_1"
    equilibrium = "rayleigh_taylor"
    geometry = "Cartesian"

    parameters = {
        "k2": 0.0,
        "k3": 1.0,
        "cte_rho0": 1.0,
        "cte_p0": 1000.0,
        "delta": -5.0,
        "g": 15.0,
        "alpha": 0.0,
        "theta": 0.35 * np.pi,
        "p1": 0.2,
        "p2": 0.6,
        "p3": 0.0,
        "p4": -0.35 * np.pi,
        "tau": 0.0,
    }
    physics_settings = {
        "flow": True,
        "external_gravity": True,
    }

    spectrum_limits = [
        {"xlim": (-1.3e4, 1.3e4), "ylim": (-1.6, 1.6)},
        {"xlim": (-500, 500), "ylim": (-1.6, 1.6)},
        {"xlim": (-10, 15), "ylim": (-1.6, 1.6)},
        {"xlim": (0.01, 0.23), "ylim": (-1.6, 1.6)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)
