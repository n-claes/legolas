from .regression import MultiRegressionTest
import numpy as np
import pytest


class TestGravito_MHD_Multi(MultiRegressionTest):
    name = "gravito MHD multi"
    filename = "multi_gravito_mhd"
    equilibrium = "gravito_mhd"
    geometry = "Cartesian"
    number_of_runs = 18
    gridpoints = 31

    kvals = np.linspace(0, np.sqrt(250), number_of_runs)

    parameters = {"k2": kvals, "k3": kvals, "cte_p0": 0.5, "g": 20.0, "alpha": 20.0}
    physics_settings = {"external_gravity": True}

    multispectrum_settings = {
        "xdata": 2 * kvals**2,  # k0**2
        "xlim": (0, 550),
        "ylim": (0, 550),
    }

    def test_gravity_value(self, series_test, series_base):
        for ds_test, ds_base in zip(series_test, series_base):
            assert np.all(
                ds_test.equilibria.get("grav") == pytest.approx(self.parameters["g"])
            )
            assert np.all(
                ds_base.equilibria.get("grav") == pytest.approx(self.parameters["g"])
            )
