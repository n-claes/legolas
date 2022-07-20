from .regression import MultiRegressionTest
import numpy as np


class TestGravitoAcousticMulti(MultiRegressionTest):
    name = "gravito acoustic multi (HD)"
    filename = "multi_gravito_acoustic"
    equilibrium = "gravito_acoustic"
    geometry = "Cartesian"
    number_of_runs = 18
    gridpoints = 31

    kvals = np.linspace(0, np.sqrt(250), number_of_runs)

    parameters = {"k2": kvals, "k3": kvals, "cte_p0": 1, "g": 0.5, "alpha": 20.42}
    physics_settings = {"external_gravity": True}

    multispectrum_settings = {
        "xdata": 2 * kvals**2,  # k0**2
        "y_scaling": 1 / 0.04080966,  # 1 / cs**2
        "xlim": (0, 550),
        "ylim": (0, 550),
    }
