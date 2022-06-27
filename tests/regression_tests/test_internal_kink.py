from .regression import RegressionTest
import pytest


class TestInternalKinkModesQR(RegressionTest):
    name = "internal kink modes k2=1 k3=0.8"
    filename = "internal_kink_QR_k2_1_k3_0.8"
    equilibrium = "internal_kink"
    geometry = "cylindrical"

    parameters = {
        "k2": 1.0,
        "k3": 0.16 * 5,
        "cte_rho0": 1.0,
        "cte_v03": 1.0,
        "cte_p0": 3.0,
        "alpha": 5.0,
    }
    physics_settings = {"flow": True}
    eigenfunction_settings = {"write_eigenfunctions": True}

    spectrum_limits = [
        {"xlim": (-7000, 7000), "ylim": (-0.08, 0.08)},
        {"xlim": (-150, 150), "ylim": (-0.08, 0.08)},
        {"xlim": (-4, 5), "ylim": (-0.08, 0.08), "RMS_TOLERANCE": 2.3},
        {"xlim": (-0.05, 1.2), "ylim": (-0.08, 0.08)},
    ]
    eigenfunctions = [
        {"eigenvalue": 0.470629 + 0.0607j},
        {"eigenvalue": 0.470629 - 0.0607j},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)
