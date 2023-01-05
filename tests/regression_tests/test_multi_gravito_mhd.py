import numpy as np
import pytest

from .regression import MultiRegressionTest


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

    def test_gravity_value(self, series_test):
        for ds_test in series_test:
            assert np.all(
                ds_test.equilibria.get("gravity") == pytest.approx(self.parameters["g"])
            )
