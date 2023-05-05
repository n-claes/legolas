import pytest

from .regression import RegressionTest

eigenfunctions = [
    {"eigenvalue": 0.20249 - 0.11790j},
    {"eigenvalue": 0.79751 - 0.11790j},
    {"eigenvalue": 0.35402 - 0.20539j},
    {"eigenvalue": 0.64598 - 0.20539j},
    {"eigenvalue": 0.48448 - 0.28639j},
    {"eigenvalue": 0.51552 - 0.28639j},
    {"eigenvalue": 0.33332 - 0.29077j},
    {"eigenvalue": 0.66668 - 0.29077j},
    {"eigenvalue": 0.5 - 0.44518624j},
]

spectrum_limits = [
    {"xlim": (-430, 430), "ylim": (-155, 5)},
    {"xlim": (-15, 15), "ylim": (-80, 10)},
    {"xlim": (-0.2, 1.2), "ylim": (-40, 10)},
    {"xlim": (-0.05, 1.1), "ylim": (-1.1, 0.2)},
]


class CouetteFlowHeating(RegressionTest):
    name = "Couette flow heating k2=0 k3=1"
    equilibrium = "couette_flow"
    geometry = "Cartesian"

    parameters = {
        "k2": 0,
        "k3": 1,
        "cte_rho0": 1.0,
        "cte_T0": 1.0,
        "cte_v02": 0.0,
        "cte_v03": 1.0,
    }
    eigenfunction_settings = {
        "write_eigenfunctions": True,
        "write_derived_eigenfunctions": True,
        "write_eigenfunction_subset": True,
        "eigenfunction_subset_center": 0.5 - 0.5j,
        "eigenfunction_subset_radius": 0.49,
    }
    physics_settings = {
        "flow": True,
        "viscosity": True,
        "viscosity_value": 1e-3,
        "viscous_heating": True,
    }


class TestCouetteFlowHeatingQR(CouetteFlowHeating):
    name = "Couette flow heating k2=0 k3=1 QR"
    filename = "couette_heating_QR_k2_0_k3_1"

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    @pytest.mark.parametrize("eigenfunction", eigenfunctions)
    def test_eigenfunction(self, eigenfunction, ds_test, ds_base):
        super().run_eigenfunction_test(eigenfunction, ds_test, ds_base)

    @pytest.mark.parametrize("derived_eigenfunction", eigenfunctions)
    def test_derived_eigenfunction(self, derived_eigenfunction, ds_test, ds_base):
        super().run_derived_eigenfunction_test(derived_eigenfunction, ds_test, ds_base)


class TestCouetteFlowHeatingQZ(CouetteFlowHeating):
    name = "Couette flow heating k2=0 k3=1 QZ"
    filename = "couette_heating_QZ_k2_0_k3_1"
    use_custom_baseline = "couette_heating_QR_k2_0_k3_1"
    solver_settings = {"solver": "QZ-direct"}

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)


class TestCouetteFlowHeatingQRCholesky(CouetteFlowHeating):
    name = "Couette flow heating k2=0 k3=1 QR Cholesky"
    filename = "couette_heating_QR_Cholesky_k2_0_k3_1"
    use_custom_baseline = "couette_heating_QR_k2_0_k3_1"
    solver_settings = {"solver": "QR-cholesky"}

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)


class TestCouetteFlowHeatingSI(CouetteFlowHeating):
    filename = "couette_heating_SI_k2_0_k3_1"
    solver_settings = {
        "solver": "arnoldi",
        "arpack_mode": "shift-invert",
        "number_of_eigenvalues": 20,
        "which_eigenvalues": "LM",
        "sigma": 0.5 - 0.6j,
    }

    spectrum_limits = [
        {"xlim": (0.1, 0.9), "ylim": (-1, 0.05)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)
