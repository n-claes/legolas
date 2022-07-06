from .regression import RegressionTest
import pytest


class RTI_ThetaPinchHD(RegressionTest):
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


class TestRTI_ThetaPinchHD_QR(RTI_ThetaPinchHD):
    name = "theta pinch HD k2=1 k3=0 QR"
    filename = "rti_theta_pinch_hd_QR_k2_1_k3_0"

    spectrum_limits = [
        {"xlim": (-1e4, 1e4), "ylim": (-0.5, 0.5)},
        {"xlim": (-100, 100), "ylim": (-0.5, 0.5)},
        {"xlim": (-0.5, 2), "ylim": (-0.5, 0.5)},
        {"xlim": (0.8, 1.08), "ylim": (-0.5, 0.5), "RMS_TOLERANCE": 3.7},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)


class TestRTI_ThetaPinchHD_SI(RTI_ThetaPinchHD):
    name = "theta pinch HD k2=1 k3=0 shift-invert"
    filename = "rti_theta_pinch_hd_SI_k2_1_k3_0"
    solver_settings = {
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 15,
        "which_eigenvalues": "LM",
        "sigma": 1.0 + 0.5j,
    }
    spectrum_limits = [
        {"xlim": (0.75, 1.15), "ylim": (-0.05, 0.5)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)


class RTI_ThetaPinchMHD(RegressionTest):
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


class TestRTI_ThetaPinch_MHD_QR(RTI_ThetaPinchMHD):
    name = "theta pinch MHD k2=1 k3=0.1 QR"
    filename = "rti_theta_pinch_mhd_QR_k2_1_k3_0.1"

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


class TestRTI_ThetaPinchMHD_SI(RTI_ThetaPinchMHD):
    name = "theta pinch MHD k2=1 k3=0 shift-invert"
    filename = "rti_theta_pinch_mhd_SI_k2_1_k3_0.1"
    solver_settings = {
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 10,
        "which_eigenvalues": "LM",
        "sigma": 1 + 0.2j,
    }
    spectrum_limits = [
        {"xlim": (0.75, 1.15), "ylim": (-0.05, 0.5)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)
