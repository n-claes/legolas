from .regression import RegressionTest
import numpy as np
import pytest


class TestUniHallAdiabaticQR(RegressionTest):
    name = "uniform adiabatic and hall k2=0.5pi k3=0.87"
    filename = "uni_adiab_hall_QR_k2_0.5pi_k3_0.87"
    equilibrium = "adiabatic_homo"
    geometry = "Cartesian"
    x_start = 0
    x_end = 1000

    parameters = {
        "k2": round(np.pi * np.sin(np.pi / 6), 14),
        "k3": round(np.pi * np.cos(np.pi / 6), 14),
        "cte_rho0": 1.0,
        "cte_T0": 1.0,
        "cte_B02": 0.0,
        "cte_B03": 1.0,
    }
    physics_settings = {
        "hall_mhd": True,
        "electron_fraction": 0.5,
        "cgs_units": True,
        "unit_density": 1.7e-14,
        "unit_magneticfield": 10,
        "unit_length": 7.534209349981049e-9,
    }
    eigenfunction_settings = {"write_eigenfunctions": True}

    eigenvalues_are_real = True
    spectrum_limits = [
        {"xlim": (-11, 11), "ylim": (-0.05, 0.05)},
        {"xlim": (0.5, 5.5), "ylim": (-0.05, 0.05)},
        {"xlim": (0.7868, 0.7878), "ylim": (-0.05, 0.05)},
        {"xlim": (4.013, 4.024), "ylim": (-0.05, 0.05)},
        {"xlim": (9.4, 9.51), "ylim": (-0.05, 0.05)},
    ]
    eigenfunctions = [
        {"eigenvalue": 0.78774536},
        {"eigenvalue": 0.78774434},
        {"eigenvalue": 0.78774264},
        {"eigenvalue": 4.0167199},
        {"eigenvalue": 9.48806988},
        {"eigenvalue": 9.48808336},
        {"eigenvalue": 9.48810584},
        {"eigenvalue": 9.4881373},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)
