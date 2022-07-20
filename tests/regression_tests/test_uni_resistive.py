from .regression import RegressionTest
import pytest
import numpy as np


class UniResistive(RegressionTest):
    equilibrium = "resistive_homo"
    geometry = "Cartesian"

    parameters = {
        "k2": 0,
        "k3": 1.0,
        "beta": 0.25,
        "cte_rho0": 1.0,
        "cte_B02": 0.0,
        "cte_B03": 1.0,
    }
    physics_settings = {
        "resistivity": True,
        "use_fixed_resistivity": True,
        "fixed_eta_value": 0.001,
    }
    eigenfunction_settings = {
        "write_eigenfunctions": True,
        "write_derived_eigenfunctions": True,
        "write_eigenfunction_subset": True,
        "eigenfunction_subset_center": 12 + 0j,
        "eigenfunction_subset_radius": 10,
    }


class TestUniResistiveQR(UniResistive):
    name = "uniform resistive k2=0 k3=1 QR"
    filename = "uni_resistive_QR_k2_0_k3_1"

    spectrum_limits = [
        {"xlim": (-375, 375), "ylim": (-110, 5)},
        {"xlim": (-30, 30), "ylim": (-30, 3)},
        {"xlim": (-10, 10), "ylim": (-7.5, 1)},
        {"xlim": (-1.2, 1.2), "ylim": (-1.1, 0.05)},
    ]
    eigenfunctions = [
        {"eigenvalue": 3.59990691 - 0.00454643j},
        {"eigenvalue": 6.98125212 - 0.01679692j},
        {"eigenvalue": 10.4098540 - 0.03721645j},
        {"eigenvalue": 13.8506386 - 0.06580399j},
        {"eigenvalue": 17.2962673 - 0.10255935j},
        {"eigenvalue": 20.7442481 - 0.14748239j},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)

    def test_eta_value(self, ds_test):
        assert np.all(
            ds_test.equilibria.get("eta")
            == pytest.approx(self.physics_settings["fixed_eta_value"])
        )


class TestUniResistiveSI(UniResistive):
    name = "uniform resistive k2=0 k3=1 shift-invert"
    filename = "uni_resistive_SI_k2_0_k3_1"
    solver_settings = {
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 30,
        "which_eigenvalues": "LM",
        "sigma": 10.0 - 0.05j,
    }

    spectrum_limits = [
        {"xlim": (-0.003, 18), "ylim": (-0.175, 0.008)},
    ]
    eigenfunctions = [
        {"eigenvalue": 3.59990691 - 0.00454643j},
        {"eigenvalue": 6.98125212 - 0.01679692j},
        {"eigenvalue": 10.4098540 - 0.03721645j},
        {"eigenvalue": 13.8506386 - 0.06580399j},
        {"eigenvalue": 17.2962673 - 0.10255935j},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)
