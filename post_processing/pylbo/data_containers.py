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
    read_ef_grid,
    read_eigenfunction,
)
from pylbo.utilities.logger import pylboLogger
from pylbo.exceptions import MatricesNotPresent
from pylbo.visualisation.continua import calculate_continua
from pylbo.utilities.toolbox import transform_to_numpy


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

    @abstractmethod
    def get_tube_speed(self, which_values=None):
        pass

    @abstractmethod
    def get_reynolds_nb(self, which_values=None):
        pass

    @abstractmethod
    def get_magnetic_reynolds_nb(self, which_values=None):
        pass

    @abstractmethod
    def get_k0_squared(self):
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

    @staticmethod
    def _get_values(array, which_values):
        if which_values is None:
            return array
        elif which_values == "average":
            return np.average(array)
        elif which_values == "minimum":
            return np.min(array)
        elif which_values == "maximum":
            return np.max(array)
        else:
            raise ValueError(f"unknown argument which_values: {which_values}")

    def get_sound_speed(self, which_values=None):
        """
        Calculates the sound speed based on the equilibrium arrays,
        given by :math:`c_s = \\sqrt{\\frac{\\gamma p_0}{\\rho_0}}`.

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
        cs : float, numpy.ndarray(dtype=float, ndim=1)
            Array with the sound speed at every grid point, or a float corresponding
            to the value of `which_values` if provided.
        """
        pressure = self.equilibria["T0"] * self.equilibria["rho0"]
        cs = np.sqrt(self.gamma * pressure / self.equilibria["rho0"])
        return self._get_values(cs, which_values)

    def get_alfven_speed(self, which_values=None):
        """
        Calculates the Alfvén speed based on the equilibrium arrays,
        given by :math:`c_A = \\sqrt{\\frac{B_0^2}{\\rho_0}}`.

        Parameters
        ----------
        which_values : str
            Which values to return:
            - None : returns the Alfvén speed as a function of the grid.
            - "average" : returns the average Alfvén speed over the grid.
            - "minimum" : returns the minimum Alfvén speed over the grid.
            - "maximum" : returns the maximum Alfvén speed over the grid

        Returns
        -------
        cA : float, numpy.ndarray(dtype=float, ndim=1)
            Array with the Alfvén speed at every grid point,
            or a float corresponding to the value of `which_values` if provided.
        """
        B0 = np.sqrt(self.equilibria["B02"] ** 2 + self.equilibria["B03"] ** 2)
        cA = B0 / np.sqrt(self.equilibria["rho0"])
        return self._get_values(cA, which_values)

    def get_tube_speed(self, which_values=None):
        """
        Calculates the tube speed for a cylinder, given by
        :math:`c_t = \\frac{c_s c_A}{\\sqrt{c_s^2 + c_A^2}}`.

        Parameters
        ----------
        which_values : str
            Which values to return:
            - None : returns the tube speed as a function of the grid.
            - "average" : returns the average tube speed over the grid.
            - "minimum" : returns the minimum tube speed over the grid.
            - "maximum" : returns the maximum tube speed over the grid

        Returns
        -------
        ct = float, numpy.ndarray(dtype=float, ndim=1)
            Array with the tube speed at every grid point,
            or a float corresponding to the value of `which_values` if provided.
            Returns `None` if the geometry is not cylindrical.
        """
        if not self.geometry == "cylindrical":
            pylboLogger.warning(
                "geometry is not cylindrical, unable to calculate tube speed"
            )
            return None
        else:
            cA = self.get_alfven_speed()
            cs = self.get_sound_speed()
            ct = cs * cA / np.sqrt(cs ** 2 + cA ** 2)
            return self._get_values(ct, which_values)

    def get_reynolds_nb(self, which_values=None):
        """
        Calculates the Reynolds number, defined as
        :math:`R_e = \\frac{ac_s}{\\eta}` where the slabsize is given by
        :math:`a = x_{end} - x_{start}`.

        Parameters
        ----------
        which_values : str
            Which values to return:
            - None : returns the Reynolds number as a function of the grid.
            - "average" : returns the average Reynolds number over the grid.
            - "minimum" : returns the minimum Reynolds number over the grid.
            - "maximum" : returns the maximum Reynolds number the grid

        Returns
        -------
        Re : float, numpy.ndarray(dtype=float, ndim=1)
            Array with the Reynolds number at every grid point,
            or a float corresponding to the value of `which_values` if provided.
            Returns `None` if the resistivity is zero somewhere on the domain.
        """
        cs = self.get_sound_speed()
        a = self.x_end - self.x_start
        eta = self.equilibria["eta"]
        if (eta == 0).any():
            pylboLogger.warning(
                "resistivity is zero somewhere on the domain, unable to "
                "calculate the Reynolds number"
            )
            return None
        else:
            Re = a * cs / eta
            return self._get_values(Re, which_values)

    def get_magnetic_reynolds_nb(self, which_values=None):
        """
        Calculates the magnetic Reynolds number, defined as
        :math:`R_m = \\frac{ac_A}{\\eta}` where the slabsize is given by
        :math:`a = x_{end} - x_{start}`.

        Parameters
        ----------
        which_values : str
            Which values to return:
            - None : returns the magnetic Reynolds number as a function of the grid.
            - "average" : returns the average magnetic Reynolds number over the grid.
            - "minimum" : returns the minimum magnetic Reynolds number over the grid.
            - "maximum" : returns the maximum magnetic Reynolds number over the grid

        Returns
        -------
        Rm : float, numpy.ndarray(dtype=float, ndim=1)
            Array with the magnetic Reynolds number at every grid point,
            or a float corresponding to the value of `which_values` if provided.
            Returns `None` if the resistivity is zero somewhere on the domain.
        """
        cA = self.get_alfven_speed()
        a = self.x_end - self.x_start
        eta = self.equilibria["eta"]
        if (eta == 0).any():
            pylboLogger.warning(
                "resistivity is zero somewhere on the domain, unable to "
                "calculate the magnetic Reynolds number"
            )
            return None
        else:
            Rm = a * cA / eta
            return self._get_values(Rm, which_values)

    def get_k0_squared(self):
        """
        Calculates the squared wave number, defined as
        :math:`k_0^2 = k_2^2 + k_3^2`.

        Returns
        -------
        k0_sq : float
            The wave number squared.

        """
        k0_sq = self.parameters.get("k2") ** 2 + self.parameters.get("k3") ** 2
        return k0_sq

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

    def get_eigenfunctions(self, ev_guesses=None, ev_idxs=None):
        """
        Returns the eigenfunctions based on given eigenvalue guesses or their
        indices. An array will be returned where every item is a dictionary, containing
        both the eigenvalue and its eigenfunctions. Either eigenvalue guesses or
        indices can be supplied, but not both.

        Parameters
        ----------
        ev_guesses : (list of) int, float, complex
            Eigenvalue guesses.
        ev_idxs : (list of) int
            Indices corresponding to the eigenvalues that need to be retrieved.

        Returns
        -------
        eigenfuncs : numpy.ndarray(dtype=dict, ndim=1)
            Array containing the eigenfunctions and eigenvalues corresponding to the
            supplied indices. Every index in this array contains a dictionary with the
            eigenfunctions and corresponding eigenvalue.
            The keys of each dictionary are the eigenfunction names.
        """
        if ev_guesses is not None and ev_idxs is not None:
            raise ValueError(
                "get_eigenfunctions: either provide guesses or indices but not both"
            )
        if ev_guesses is not None:
            idxs, _ = self.get_nearest_eigenvalues(ev_guesses)
        else:
            idxs = transform_to_numpy(ev_idxs)
            for idx in idxs:
                if not isinstance(idx, (int, np.int64)):
                    raise ValueError("get_eigenfunctions: ev_idxs should be integers")
        eigenfuncs = np.array([{}] * len(idxs), dtype=dict)
        with open(self.datfile, "rb") as istream:
            for dict_idx, ef_idx in enumerate(idxs):
                efs = read_eigenfunction(istream, self.header, ef_idx)
                efs.update({"eigenvalue": self.eigenvalues[ef_idx]})
                eigenfuncs[dict_idx] = efs
        return eigenfuncs

    def get_nearest_eigenvalues(self, ev_guesses):
        """
        Calculates the eigenvalues nearest to a given guess. This calculates
        the nearest eigenvalue based on the distance between two points.

        Parameters
        ----------
        ev_guesses : float, complex, list of float, list of complex
            The guesses for the eigenvalues. These can be a single float/complex value,
            or a list/Numpy array of floats/complex values.

        Returns
        -------
        idxs : numpy.ndarray(dtype=int, ndim=1)
            The indices of the nearest eigenvalues in the :attr:`eigenvalues` array.
        eigenvalues: numpy.ndarray(dtype=complex, ndim=1)
            The nearest eigenvalues to the provided guesses, corresponding with the
            indices `idxs`.
        """
        ev_guesses = transform_to_numpy(ev_guesses)
        idxs = np.empty(shape=len(ev_guesses), dtype=np.int)
        eigenvals = np.empty(shape=len(ev_guesses), dtype=np.complex)
        for i, ev_guess in enumerate(ev_guesses):
            # distance from guess to all eigenvalues
            distances = (self.eigenvalues.real - ev_guess.real) ** 2 + (
                self.eigenvalues.imag - ev_guess.imag
            ) ** 2
            # closest distance (squared)
            idx = distances.argmin()
            idxs[i] = idx
            eigenvals[i] = self.eigenvalues[idx]
        return idxs, eigenvals


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
            return self.datasets[idx.start : idx.stop : idx.step]
        else:
            return self.datasets[idx]

    def __len__(self):
        return len(self.datasets)

    @property
    def efs_written(self):
        return np.array([ds.efs_written for ds in self.datasets])

    @property
    def ef_names(self):
        names = np.array([ds.ef_names for ds in self.datasets], dtype=object)
        for item in names:
            if item is None:
                continue
            else:
                return item
        else:
            return None

    @property
    def ef_grid(self):
        return np.array([ds.ef_grid for ds in self.datasets], dtype=object)

    def get_sound_speed(self, which_values=None):
        return np.array([ds.get_sound_speed(which_values) for ds in self.datasets])

    def get_alfven_speed(self, which_values=None):
        return np.array([ds.get_alfven_speed(which_values) for ds in self.datasets])

    def get_tube_speed(self, which_values=None):
        return np.array([ds.get_tube_speed(which_values) for ds in self.datasets])

    def get_reynolds_nb(self, which_values=None):
        return np.array([ds.get_reynolds_nb(which_values) for ds in self.datasets])

    def get_magnetic_reynolds_nb(self, which_values=None):
        return np.array(
            [ds.get_magnetic_reynolds_nb(which_values) for ds in self.datasets]
        )

    def get_k0_squared(self):
        return np.array([ds.get_k0_squared() for ds in self.datasets])
