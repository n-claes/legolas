import pytest

from .regression import RegressionTest

spectrum_limits = [
    {"xlim": (-45, 45), "ylim": (-0.7, 0.7)},
    {"xlim": (-10, 10), "ylim": (-0.7, 0.7)},
    {"xlim": (-1.5, 1.5), "ylim": (-0.7, 0.7)},
    {"xlim": (-0.05, 0.05), "ylim": (-0.7, 0.7)},
    {"xlim": (-0.01, 0.01), "ylim": (-0.01, 0.7)},
]


class MRI_Accretion(RegressionTest):
    equilibrium = "MRI_accretion"
    geometry = "cylindrical"
    x_start = 1.0
    x_end = 2.0

    parameters = {"k2": 0.0, "k3": 70.0, "beta": 100.0, "tau": 1.0, "nu": 0.1}
    physics_settings = {"flow": True, "external_gravity": True}
    eigenfunction_settings = {
        "write_eigenfunctions": True,
        "write_derived_eigenfunctions": True,
        "write_eigenfunction_subset": True,
        "eigenfunction_subset_center": -0.001 + 0.59j,
        "eigenfunction_subset_radius": 0.06,
    }


class TestMRI_AccretionQR(MRI_Accretion):
    name = "MRI accretion k2=0 k3=70 QZ"
    filename = "mri_accretion_QR_k2_0_k3_70"

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

    @pytest.mark.parametrize("derived_eigenfunction", eigenfunctions)
    def test_derived_eigenfunction(self, derived_eigenfunction, ds_test, ds_base):
        super().run_derived_eigenfunction_test(derived_eigenfunction, ds_test, ds_base)


class TestMRI_AccretionQZ(MRI_Accretion):
    name = "MRI accretion k2=0 k3=70 QZ"
    filename = "mri_accretion_QZ_k2_0_k3_70"
    use_custom_baseline = "mri_accretion_QR_k2_0_k3_70"
    solver_settings = {"solver": "QZ-direct"}

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)


class TestMRI_AccretionQRCholesky(MRI_Accretion):
    name = "MRI accretion k2=0 k3=70 QR Cholesky"
    filename = "mri_accretion_QR_cholesky_k2_0_k3_70"
    use_custom_baseline = "mri_accretion_QR_k2_0_k3_70"
    solver_settings = {"solver": "QR-cholesky"}

    eigenfunctions = [
        {"eigenvalue": -0.002028 + 0.627722j, "RMS_TOLERANCE": 3.8},
        {"eigenvalue": -0.00186038 + 0.5804956j, "RMS_TOLERANCE": 4.7},
        {"eigenvalue": -0.0017397 + 0.5443967j, "RMS_TOLERANCE": 5.6},
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
