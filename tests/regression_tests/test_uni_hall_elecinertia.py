import numpy as np
import pytest

from .regression import RegressionTest


class UniHallElectronInertia(RegressionTest):
    equilibrium = "adiabatic_homo"
    geometry = "Cartesian"
    x_start = 0
    x_end = 1000

    parameters = {
        "k2": 10 * np.sin(np.pi / 6),
        "k3": 10 * np.cos(np.pi / 6),
        "cte_rho0": 1,
        "cte_T0": 1,
        "cte_B02": 0,
        "cte_B03": 1,
    }
    physics_settings = {
        "hall_mhd": True,
        "hall_substitution": True,
        "electron_fraction": 0.5,
        "elec_inertia": True,
        "unit_density": 1.7e-14,
        "unit_magneticfield": 10,
        "unit_length": 7.534209349981049e-9,
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
        inertia = 0.00054462
        assert np.allclose(ds_test.equilibria["inertia"], inertia)
        assert np.allclose(ds_base.equilibria["inertia"], inertia)


class TestUniHallElectronEnertiaQR(UniHallElectronInertia):
    name = "uniform hall k2=5 k3=8.66 QR"
    filename = "uni_hall_elecinertia_QR_k2_5_k3_8.66"
    spectrum_limits = [
        {"xlim": (-85, 85), "ylim": (-0.05, 0.05)},
        {"xlim": (-15, 15), "ylim": (-0.05, 0.05)},
        {"xlim": (79.5, 83.5), "ylim": (-0.05, 0.05)},
        {"xlim": (10, 15), "ylim": (-0.05, 0.05)},
        {"xlim": (0.85, 1), "ylim": (-0.05, 0.05)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)


class TestUniHallElectronInertiaSI(UniHallElectronInertia):
    name = "uniform hall k2=5 k3=8.66 shift-invert"
    filename = "uni_hall_elecinertia_SI_k2_5_k3_8.66"
    solver_settings = {
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 10,
        "which_eigenvalues": "LM",
        "sigma": 9.44805 + 0j,
    }
    spectrum_limits = [
        {"xlim": (7, 13), "ylim": (-0.05, 0.05)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)


class TestUniHallElectronInertiaSI_EFS(UniHallElectronInertia):
    name = "uniform hall k2=5 k3=8.66 shift-invert efs"
    filename = "uni_hall_elecinertia_SI2_k2_5_k3_8.66"
    solver_settings = {
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 10,
        "which_eigenvalues": "LM",
        "sigma": 83.1316 + 0j,
    }
    spectrum_limits = [
        {"xlim": (83.131, 83.132), "ylim": (-0.05, 0.05)},
    ]
    eigenfunction_settings = {
        "write_eigenfunctions": True,
    }
    eigenfunctions = [
        {"eigenvalue": 83.13145786},
        {"eigenvalue": 83.13146882},
        {"eigenvalue": 83.13148708},
        {"eigenvalue": 83.13151265},
        {"eigenvalue": 83.13154552},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)
