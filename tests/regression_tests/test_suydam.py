from .regression import RegressionTest
import pytest


class TestSuydamQR(RegressionTest):
    name = "Suydam modes k2=1 k3=-1.2"
    filename = "suydam_QR_k2_1_k3_-1.2"
    equilibrium = "suydam_cluster"
    geometry = "cylindrical"

    parameters = {
        "k2": 1.0,
        "k3": -1.2,
        "cte_rho0": 1,
        "cte_v02": 0,
        "cte_v03": 0.14,
        "cte_p0": 0.05,
        "p1": 0.1,
        "alpha": 2.0,
    }

    physics_settings = {"flow": True}

    spectrum_limits = [
        {"xlim": (-350, 350), "ylim": (-0.01, 0.01), "RMS_TOLERANCE": 3},
        {"xlim": (-1, 1), "ylim": (-0.01, 0.01)},
        {"xlim": (-0.2, -0.02), "ylim": (-0.01, 0.01), "RMS_TOLERANCE": 3.15},
        {"xlim": (-0.16, -0.08), "ylim": (-0.01, 0.01), "RMS_TOLERANCE": 2.2},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)
