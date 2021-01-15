import numpy as np
from pathlib import Path
from abc import ABC, abstractmethod
from pylbo.utilities.datfile_utils import (
    get_header,
    read_grid,
    read_grid_gauss,
    read_eigenvalues,
    read_equilibrium_arrays,
    read_matrix_B,
    read_matrix_A,
)
from pylbo.utilities.logger import pylboLogger
from pylbo.exceptions import MatricesNotPresent
from pylbo.continua import calculate_continua


class LegolasDataContainer(ABC):
    @property
    def continua(self):
        if isinstance(self, LegolasDataSet):
            return self._continua
        elif isinstance(self, LegolasDataSeries):
            keys = self[0].continua.keys()
            _continua = {key: [] for key in keys}
            for ds in self:
                for key in keys:
                    _continua[key].append(ds.continua[key])
            return {key: np.array(values) for key, values in _continua.items()}
        else:
            raise TypeError(f"unexpected instance: {type(self)}")

    @property
    def parameters(self):
        if isinstance(self, LegolasDataSet):
            return self._parameters
        elif isinstance(self, LegolasDataSeries):
            keys = self[0].parameters.keys()
            _params = {key: [] for key in keys}
            for ds in self:
                for key in keys:
                    _params[key].append(ds.parameters[key])
            return {key: np.array(values) for key, values in _params.items()}
        else:
            raise TypeError(f"unexpected instance: {type(self)}")

    @abstractmethod
    def get_sound_speed(self, **kwargs):
        pass


class LegolasDataSet(LegolasDataContainer):
    def __init__(self, datfile):
        self.datfile = Path(datfile)
        with open(self.datfile, "rb") as istream:
            self.header = get_header(istream)
            self.grid = read_grid(istream, self.header)
            self.grid_gauss = read_grid_gauss(istream, self.header)
            self.equilibria = read_equilibrium_arrays(istream, self.header)
            self.eigenvalues = read_eigenvalues(istream, self.header)

        self.geometry = self.header["geometry"]
        if self.geometry == "Cartesian":
            self.scale_factor = np.ones_like(self.grid_gauss)
            pylboLogger.debug("dataset: scale factor set to unity.")
        else:
            self.scale_factor = self.grid_gauss
            pylboLogger.debug("dataset: scale factor set to radial coordinate.")

        self.x_start = self.header["x_start"]
        self.x_end = self.header["x_end"]
        self.gridpoints = self.header["gridpts"]
        self.gauss_gridpoints = self.header["gauss_gridpts"]
        self.matrix_gridpoints = self.header["matrix_gridpts"]
        self.ef_gridpoints = self.header["ef_gridpts"]
        self.gamma = self.header["gamma"]
        self.eq_type = self.header["eq_type"]
        self._parameters = self.header["params"]
        self.cgs = self.header["cgs"]
        self.units = self.header["units"]
        self.eq_names = self.header["equil_names"]
        self.legolas_version = self.header["legolas_version"]

        self._continua = calculate_continua(self)

    def __iter__(self):
        yield self

    def get_sound_speed(self):
        """
        Calculates the sound speed based on the equilibrium arrays.

        Returns
        -------
        cs : numpy.ndarray(dtype=float, ndim=1)
            Sound speed at every grid point, defined as
            :math:`c_s = \\sqrt{\\frac{\\gamma p_0}{\\rho_0}}`.
        """
        pressure = self.equilibria["T0"] * self.equilibria["rho0"]
        cs = np.sqrt(self.gamma * pressure / self.equilibria["rho0"])
        return cs

    def get_matrix_B(self):
        """
        Retrieves the matrix B from the datfile.

        Returns
        -------
        3-tuple of numpy.ndarray(dtype=int, ndim=1)
            Tuple containing the rows, columns and values of the non-zero B-matrix
            elements. Rows and columns are integers, values are real.

        Raises
        ------
        MatricesNotPresent
            If the matrices were not saved to the datfile.
        """
        if not self.header["matrices_written"]:
            raise MatricesNotPresent(self.datfile)
        with open(self.datfile, "rb") as istream:
            rows, cols, vals = read_matrix_B(istream, self.header)
        return rows, cols, vals

    def get_matrix_A(self):
        """
        Retrieves the matrix A from the datfile.

        Returns
        -------
        3-tuple of numpy.ndarray(dtype=int, ndim=1)
            Tuple containing the rows, columns and values of the non-zero A-matrix
            elements. Rows and columns are integers, values are complex.

        Raises
        ------
        MatricesNotPresent
            If the matrices were not saved to the datfile.
        """
        if not self.header["matrices_written"]:
            raise MatricesNotPresent(self.datfile)
        with open(self.datfile, "rb") as istream:
            rows, cols, vals = read_matrix_A(istream, self.header)
        return rows, cols, vals


class LegolasDataSeries(LegolasDataContainer):
    def __init__(self, datfiles):
        self.datasets = [LegolasDataSet(datfile) for datfile in datfiles]

    def __iter__(self):
        for ds in self.datasets:
            yield ds

    def __getitem__(self, idx):
        return self.datasets[idx]

    def __len__(self):
        return len(self.datasets)

    def get_sound_speed(self, which_values="average"):
        pass
