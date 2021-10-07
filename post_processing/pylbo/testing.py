import pickle
import numpy as np
from pylbo.data_containers import LegolasDataSet


def pickle_dataseries_to_file(series, filepath):
    with open(filepath, "wb") as ostream:
        pickle.dump(series, ostream, pickle.HIGHEST_PROTOCOL)


def load_pickled_dataseries(filepath):
    with open(filepath, "rb") as istream:
        series = pickle.load(istream)
    return series


class FakeDataSet(LegolasDataSet):
    def __init__(self, datfile, seed=None):
        if seed is None:
            seed = 12345
        np.random.seed(seed)
        super().__init__(datfile)

        self._set_test_equilibria()
        self._set_test_eigenvalues()

    def _set_test_equilibria(self):
        for name in self.eq_names:
            x0 = np.random.randint(0, 10)
            x1 = np.random.randint(10, 25)
            self.equilibria[name] = np.linspace(x0, x1, self.gauss_gridpoints)
            if np.random.choice([True, False]):
                self.equilibria[name] = np.flip(self.equilibria[name])

    def _set_test_eigenvalues(self):
        self.eigenvalues = (
            np.random.randn(len(self.eigenvalues), 2).view(complex).flatten()
        )
