import numpy as np

from .regression import MultiRegressionTest


class TestGravitoAcousticMultiHD(MultiRegressionTest):
    name = "gravito acoustic multi (pure HD)"
    filename = "multi_gravito_acoustic_HD"
    equilibrium = "gravito_acoustic"
    geometry = "Cartesian"
    number_of_runs = 18
    gridpoints = 31

    kvals = np.linspace(0, np.sqrt(250), number_of_runs)

    parameters = {"k2": kvals, "k3": kvals, "cte_p0": 1, "g": 0.5, "alpha": 20.42}
    physics_settings = {"physics_type": "hd", "external_gravity": True}

    multispectrum_settings = {
        "xdata": 2 * kvals**2,  # k0**2
        "y_scaling": 1 / 0.04080966,  # 1 / cs**2
        "xlim": (0, 550),
        "ylim": (0, 550),
    }

    def test_physics_type(self, series_test):
        for ds_test in series_test:
            assert ds_test.header.get("physics_type") == "hd"
