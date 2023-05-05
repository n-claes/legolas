from .regression import RegressionTest
import pytest


class TestTokamakQR(RegressionTest):
    name = "tokamak constant current k2=-2 k3=0.2"
    filename = "tokamak_cte_current_QR_k2_-2_k3_0.2"
    equilibrium = "constant_current_tokamak"
    geometry = "cylindrical"

    parameters = {
        "k2": -2.0,
        "k3": 0.2,
        "j0": (2.0 * 0.2) / 1.95,
        "cte_rho0": 1.0,
        "cte_B03": 1.0,
    }
    eigenfunction_settings = {"write_eigenfunctions": True}

    spectrum_limits = [
        {"xlim": (-350, 350), "ylim": (-0.01, 0.01), "RMS_TOLERANCE": 2.2},
        {"xlim": (-0.01, 0.01), "ylim": (-0.01, 0.01), "RMS_TOLERANCE": 3.81},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)
