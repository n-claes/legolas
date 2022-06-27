from .regression import RegressionTest
import numpy as np
import pytest


class TestUniAdiabaticQR(RegressionTest):
    name = "uniform adiabatic k2=0 k3=pi"
    filename = "uni_adiab_QR_k2_0_k3_pi"
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
        "write_eigenfunction_subset": True,
        "eigenfunction_subset_center": 20 + 1j,
        "eigenfunction_subset_radius": 15,
    }

    eigenvalues_are_real = True
    spectrum_limits = [
        {"xlim": (-600, 600), "ylim": (-0.05, 0.05)},
        {"xlim": (-50, 50), "ylim": (-0.05, 0.05)},
        {"xlim": (-0.5, 5), "ylim": (-0.05, 0.05)},
    ]

    eigenfunctions = [
        {"eigenvalue": 11.18509},
        {"eigenvalue": 16.02714},
        {"eigenvalue": 21.00395},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)
