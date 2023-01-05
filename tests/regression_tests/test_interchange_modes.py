import numpy as np
import pytest

from .regression import RegressionTest


class TestInterchangeModesQR(RegressionTest):
    name = "interchange modes k2=pi k3=pi"
    filename = "interchange_modes_QR_k2_pi_k3_pi"
    equilibrium = "interchange_modes"
    geometry = "Cartesian"

    parameters = {
        "k2": np.pi,
        "k3": np.pi,
        "g": 0.5,
        "cte_p0": 0.25,
        "lambda": 0,
        "alpha": 20.0,
    }
    physics_settings = {"external_gravity": True}
    spectrum_limits = [
        {"xlim": (-75, 75), "ylim": (-0.6, 0.6)},
        {"xlim": (-5, 5), "ylim": (-0.6, 0.6)},
        {"xlim": (-1, 1), "ylim": (-0.6, 0.6)},
        {"xlim": (-0.5, 0.5), "ylim": (-0.6, 0.6)},
    ]

    @pytest.mark.parametrize("limits", spectrum_limits)
    def test_spectrum(self, limits, ds_test, ds_base):
        super().run_spectrum_test(limits, ds_test, ds_base)

    def test_external_gravity(self, ds_test):
        assert np.all(
            ds_test.equilibria.get("gravity") == pytest.approx(self.parameters["g"])
        )
