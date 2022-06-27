from .regression import RegressionTest
import pytest


class TestFluxtubeCoronalQR(RegressionTest):
    name = "fluxtube coronal k2=0 k3=4"
    filename = "fluxtube_coronal_QR_k2_0_k3_4"
    equilibrium = "coronal_flux_tube"
    geometry = "cylindrical"
    x_start = 0
    x_end = 10

    parameters = {"k2": 0, "k3": 4.0, "cte_rho0": 1.0, "cte_p0": 1.0, "r0": 1.0}
    eigenvalues_are_real = True
    spectrum_limits = [
        {"xlim": (-6300, 6300), "ylim": (-0.05, 0.05)},
        {"xlim": (-1000, 1000), "ylim": (-0.05, 0.05), "RMS_TOLERANCE": 2.7},
        {"xlim": (-50, 50), "ylim": (-0.05, 0.05), "RMS_TOLERANCE": 2.7},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)


class TestFluxtubePhotosphericQR(RegressionTest):
    name = "fluxtube photospheric k2=0 k3=4"
    filename = "fluxtube_photospheric_QR_k2_0_k3_4"
    equilibrium = "photospheric_flux_tube"
    geometry = "cylindrical"
    x_start = 0
    x_end = 10

    parameters = {"k2": 0, "k3": 4.0, "cte_rho0": 1.0, "cte_p0": 1.0, "r0": 1.0}
    eigenvalues_are_real = True
    spectrum_limits = [
        {"xlim": (-1800, 1800), "ylim": (-0.05, 0.05)},
        {"xlim": (-600, 600), "ylim": (-0.05, 0.05), "RMS_TOLERANCE": 2.7},
        {"xlim": (-15, 15), "ylim": (-0.05, 0.05), "RMS_TOLERANCE": 2.7},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)
