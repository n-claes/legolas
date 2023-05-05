import pytest

from .regression import RegressionTest


class TestKelvinHelmholtzHD_QR(RegressionTest):
    name = "kelvin helmholtz HD k2=0 k3=1 QR"
    filename = "kelvin_helmholtz_HD_QR_k2_0_k3_1"
    equilibrium = "kelvin_helmholtz"
    geometry = "Cartesian"

    parameters = {
        "k2": 0.0,
        "k3": 1.0,
        "cte_rho0": 1.0,
        "cte_p0": 10.0,
        "delta": 0.0,
        "g": 0.0,
        "alpha": 0.0,
        "theta": 0.0,
        "p1": 0.0,
        "p2": 0.0,
        "p3": 1.0,
        "p4": 0.0,
        "tau": 11.0,
    }
    physics_settings = {"physics_type": "hd", "flow": True}

    spectrum_limits = [
        {"xlim": (-300, 300), "ylim": (-0.2, 0.2)},
        {"xlim": (-30, 30), "ylim": (-0.2, 0.2)},
        {"xlim": (-1.2, 1.2), "ylim": (-0.2, 0.2)},
    ]

    @pytest.mark.required
    def test_physics_type(self, ds_test):
        assert ds_test.header.get("physics_type") == "hd"

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)
