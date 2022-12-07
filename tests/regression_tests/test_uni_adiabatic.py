import numpy as np
import pytest

from .regression import RegressionTest

eigenfunctions = [
    {"eigenvalue": 6.745518},
    {"eigenvalue": 11.18509},
    {"eigenvalue": 16.02714},
    {"eigenvalue": 21.00395},
]


class UniAdiabatic(RegressionTest):
    equilibrium = "adiabatic_homo"
    geometry = "Cartesian"

    parameters = {
        "k2": 0,
        "k3": np.pi,
        "cte_rho0": 1.0,
        "cte_T0": 1.0,
        "cte_B02": 0.0,
        "cte_B03": 1.0,
    }
    eigenfunction_settings = {
        "write_eigenfunctions": True,
        "write_derived_eigenfunctions": True,
        "write_eigenfunction_subset": True,
        "eigenfunction_subset_center": 20 + 1j,
        "eigenfunction_subset_radius": 15,
    }
    eigenvalues_are_real = True


class TestUniAdiabaticQR(UniAdiabatic):
    name = "uniform adiabatic k2=0 k3=pi QR"
    filename = "uni_adiab_QR_k2_0_k3_pi"

    spectrum_limits = [
        {"xlim": (-600, 600), "ylim": (-0.05, 0.05)},
        {"xlim": (-50, 50), "ylim": (-0.05, 0.05)},
        {"xlim": (-0.5, 5), "ylim": (-0.05, 0.05)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)

    @pytest.mark.parametrize("derived_eigenfunction", eigenfunctions)
    def test_derived_eigenfunction(self, derived_eigenfunction, ds_test, ds_base):
        super().run_derived_eigenfunction_test(derived_eigenfunction, ds_test, ds_base)


class TestUniAdiabaticQZ(TestUniAdiabaticQR):
    name = "uniform adiabatic k2=0 k3=pi QZ"
    filename = "uni_adiab_QZ_k2_0_k3_pi"
    use_custom_baseline = "uni_adiab_QR_k2_0_k3_pi"
    solver_settings = {"solver": "QZ-direct"}
    custom_evs_all_real_tol = 5e-9


class TestUniAdiabaticQRCholesky(TestUniAdiabaticQR):
    name = "uniform adiabatic k2=0 k3=pi QR Cholesky"
    filename = "uni_adiab_QR_cholesky_k2_0_k3_pi"
    use_custom_baseline = "uni_adiab_QR_k2_0_k3_pi"
    solver_settings = {"solver": "QR-cholesky"}


class TestUniAdiabaticSI(UniAdiabatic):
    name = "uniform adiabatic k2=0 k3=pi shift-invert"
    filename = "uni_adiab_SI_k2_0_k3_pi"
    solver_settings = {
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 6,
        "which_eigenvalues": "LM",
        "sigma": 15 + 0j,
    }

    spectrum_limits = [
        {"xlim": (-0.1, 30), "ylim": (-1e-5, 1e-5)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)
