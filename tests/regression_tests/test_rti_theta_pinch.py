from .regression import RegressionTest
import pytest


class Test_RTI_ThetaPinch_HD_QR(RegressionTest):
    name = "theta pinch HD k2=1 k3=0"
    filename = "rti_theta_pinch_hd_QR_k2_1_k3_0"
    equilibrium = "RTI_theta_pinch"
    geometry = "cylindrical"

    parameters = {
        "k2": 1.0,
        "k3": 0,
        "cte_rho0": 1.0,
        "alpha": 2.0,
        "delta": 1 / 6,
        "r0": 0.0,
    }

    physics_settings = {"flow": True}

    spectrum_limits = [
        {"xlim": (-1e4, 1e4), "ylim": (-0.5, 0.5)},
        {"xlim": (-100, 100), "ylim": (-0.5, 0.5)},
        {"xlim": (-0.5, 2), "ylim": (-0.5, 0.5)},
        {"xlim": (0.8, 1.08), "ylim": (-0.5, 0.5), "RMS_TOLERANCE": 3.7},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)


class Test_RTI_ThetaPinch_MHD_QR(RegressionTest):
    name = "theta pinch MHD k2=1 k3=0.1"
    filename = "rti_theta_pinch_mhd_QR_k2_1_k3_0.1"
    equilibrium = "RTI_theta_pinch"
    geometry = "cylindrical"

    parameters = {
        "k2": 1.0,
        "k3": 0.1,
        "cte_rho0": 1.0,
        "alpha": 2.0,
        "delta": 1 / 6,
        "r0": 0.0,
    }

    physics_settings = {"flow": True}

    spectrum_limits = [
        {"xlim": (-1e4, 1e4), "ylim": (-0.35, 0.35)},
        {"xlim": (-200, 200), "ylim": (-0.35, 0.35)},
        {"xlim": (-2.5, 5), "ylim": (-0.35, 0.35)},
        {"xlim": (0.05, 1.2), "ylim": (-0.35, 0.35)},
        {"xlim": (0.85, 1.06), "ylim": (-0.35, 0.35), "RMS_TOLERANCE": 2.9},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)
