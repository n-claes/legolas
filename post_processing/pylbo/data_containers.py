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
    read_ef_grid
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
    def efs_written(self):
        pass

    @abstractmethod
    def ef_grid(self):
        pass

    @abstractmethod
    def ef_names(self):
        pass

    @abstractmethod
    def get_sound_speed(self, which_values=None):
        pass

    @abstractmethod
    def get_alfven_speed(self, which_values=None):
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

    @property
    def efs_written(self):
        return self.header["eigenfuncs_written"]

    @property
    def ef_names(self):
        return self.header.get("ef_names", None)

    @property
    def ef_grid(self):
        if self.efs_written:
            with open(self.datfile, "rb") as istream:
                grid = read_ef_grid(istream, self.header)
            return grid
        else:
            return None

    def get_sound_speed(self, which_values=None):
        """
        Calculates the sound speed based on the equilibrium arrays.

        Parameters
        ----------
        which_values : str
            Which values to return:
            - None : returns the sound speed as a function of the grid.
            - "average" : returns the average sound speed over the grid.
            - "minimum" : returns the minimum sound speed over the grid.
            - "maximum" : returns the maximum sound speed over the grid

        Returns
        -------
        cs : numpy.ndarray(dtype=float, ndim=1)
            Sound speed at every grid point, defined as
            :math:`c_s = \\sqrt{\\frac{\\gamma p_0}{\\rho_0}}`.
        """
        pressure = self.equilibria["T0"] * self.equilibria["rho0"]
        cs = np.sqrt(self.gamma * pressure / self.equilibria["rho0"])
        if which_values == "average":
            return np.average(cs)
        elif which_values == "minimum":
            return np.min(cs)
        elif which_values == "maximum":
            return np.max(cs)
        return cs

    def get_alfven_speed(self, which_values=None):
        """
        Calculates the Alfvén speed based on the equilibrium arrays,
        given by :math:`c_A = \\sqrt{\\frac{B_0^2}{\\rho_0}}`.

        Parameters
        ----------
        which_values : str
            Which values to return:
            - None : returns the sound speed as a function of the grid.
            - "average" : returns the average sound speed over the grid.
            - "minimum" : returns the minimum sound speed over the grid.
            - "maximum" : returns the maximum sound speed over the grid

        Returns
        -------
        cA : numpy.ndarray(dtype=float, ndim=1)
            The Alfvén speed at every gridpoint.
        """
        B0 = np.sqrt(self.equilibria["B02"] ** 2 + self.equilibria["B03"] ** 2)
        cA = B0 / np.sqrt(self.equilibria["rho0"])
        if which_values == "average":
            return np.average(cA)
        elif which_values == "minimum":
            return np.min(cA)
        elif which_values == "maximum":
            return np.max(cA)
        return cA

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
        self.geometry = set([ds.geometry for ds in self.datasets])
        if len(self.geometry) == 1:
            self.geometry = self.geometry.pop()

    def __iter__(self):
        for ds in self.datasets:
            yield ds

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return self.datasets[idx.start:idx.stop:idx.step]
        else:
            return self.datasets[idx]

    def __len__(self):
        return len(self.datasets)

    @property
    def efs_written(self):
        return np.array([ds.efs_written for ds in self.datasets])

    @property
    def ef_names(self):
        names = np.array([ds.ef_names for ds in self.datasets])
        for item in names:
            if item is None:
                continue
            else:
                return item
        else:
            return None

    @property
    def ef_grid(self):
        return np.array([ds.ef_grid for ds in self.datasets])

    def get_sound_speed(self, which_values=None):
        return np.array([ds.get_sound_speed(which_values) for ds in self.datasets])

    def get_alfven_speed(self, which_values=None):
        return np.array([ds.get_alfven_speed(which_values) for ds in self.datasets])
