import numpy as np
import pytest

from .regression import RegressionTest


class UniHallAdiabatic(RegressionTest):
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
        "hall_substitution": True,
        "electron_fraction": 0.5,
        "unit_density": 1.7e-14,
        "unit_magneticfield": 10,
        "unit_length": 7.534209349981049e-9,
    }
    eigenfunction_settings = {
        "write_eigenfunctions": True,
        "write_eigenfunction_subset": True,
        "eigenfunction_subset_center": 5.13794 + 0j,
        "eigenfunction_subset_radius": 4.351,
    }
    eigenvalues_are_real = True

    def test_units(self, ds_test):
        assert ds_test.cgs
        for val in ("density", "magneticfield", "length"):
            assert (
                ds_test.units.get(f"unit_{val}") == self.physics_settings[f"unit_{val}"]
            )

    def test_electron_fraction(self, ds_test, ds_base):
        key = "electronfraction"
        assert (
            ds_test.parameters[key]
            == pytest.approx(ds_base.parameters[key])
            == pytest.approx(0.5)
        )

    def test_electron_inertia(self, ds_test, ds_base):
        assert np.allclose(ds_test.equilibria["inertia"], 0)
        assert np.allclose(ds_base.equilibria["inertia"], 0)


class TestUniHallAdiabaticQR(UniHallAdiabatic):
    name = "uniform adiabatic and hall k2=0.5pi k3=0.87 QR"
    filename = "uni_adiab_hall_QR_k2_0.5pi_k3_0.87"
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


class TestUniHallAdiabaticSI(UniHallAdiabatic):
    name = "uniform adiabatic and hall k2=0.5pi k3=0.87 shift-invert"
    filename = "uni_adiab_hall_SI_k2_0.5pi_k3_0.87"
    solver_settings = {
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 25,
        "which_eigenvalues": "LM",
        "sigma": 9.44805 + 0j,
    }

    eigenfunctions = [
        {"eigenvalue": 9.48806988},
        {"eigenvalue": 9.48808336},
        {"eigenvalue": 9.48810584},
        {"eigenvalue": 9.4881373},
        {"eigenvalue": 9.48817775},
        {"eigenvalue": 9.4882272},
    ]
    spectrum_limits = [
        {"xlim": (9.4875, 9.4895), "ylim": (-0.05, 0.05)},
        {"xlim": (9.4, 9.5), "ylim": (-0.05, 0.05)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)
