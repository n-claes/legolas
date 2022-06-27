from .regression import RegressionTest
import pytest


class TestMRI_AccretionQR(RegressionTest):
    name = "MRI accretion k2=0 k3=70"
    filename = "mri_accretion_QR_k2_0_k3_70"
    equilibrium = "MRI_accretion"
    geometry = "cylindrical"
    x_start = 1.0
    x_end = 2.0

    parameters = {"k2": 0.0, "k3": 70.0, "beta": 100.0, "tau": 1.0, "nu": 0.1}
    physics_settings = {"flow": True, "external_gravity": True}
    eigenfunction_settings = {"write_eigenfunctions": True}

    spectrum_limits = [
        {"xlim": (-45, 45), "ylim": (-0.7, 0.7)},
        {"xlim": (-10, 10), "ylim": (-0.7, 0.7)},
        {"xlim": (-1.5, 1.5), "ylim": (-0.7, 0.7)},
        {"xlim": (-0.05, 0.05), "ylim": (-0.7, 0.7)},
        {"xlim": (-0.01, 0.01), "ylim": (-0.01, 0.7)},
    ]
    eigenfunctions = [
        {"eigenvalue": -0.002028 + 0.627722j},
        {"eigenvalue": -0.00186038 + 0.5804956j},
        {"eigenvalue": -0.0017397 + 0.5443967j},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)
