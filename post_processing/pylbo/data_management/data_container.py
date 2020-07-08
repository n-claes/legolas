import numpy as np
from pathlib import Path
from ..utilities.datfile_utils import \
    get_header, \
    read_grid, \
    read_grid_gauss, \
    read_eigenvalues, \
    read_equilibrium_arrays, \
    read_matrix_B, \
    read_matrix_A
from ..utilities.exceptions import \
    InvalidLegolasFile, \
    MatricesNotPresent
from ..utilities.continua import get_continuum_regions
from ..utilities.logger import pylboLogger


class LegolasDataContainer:
    """
    Class acting as a general data container.

    Parameters
    ----------
    datfile : str or ~os.PathLike
        Path to a datfile from the Legolas code.

    Attributes
    ----------
    datfile : str
        Path to the datfile.
    geometry : str
        The geometry.
    x_start : float
        Start of the grid.
    x_end : float
        End of the grid.
    gridpts : int
        Amount of gridpoints in the base grid.
    gauss_gridpts : int
        Amount of gridpoints in the Gaussian grid.
    matrix_gridpts : int
        Size of the matrices A and B.
    ef_gridpts : int
        Size of an eigenfunction array.
    gamma : float
        Ratio of specific heats.
    eq_type : str
        Chosen equilibrium in the datfile.
    parameters : dict
        Various parameters such as `k2`, `k3`, `cte_rho`, etc.
    cgs : bool
        Unit system in cgs or not.
    units : dict
        Dictionary containing the unit normalisations.
    eq_names : numpy.ndarray(dtype=str, ndim=1)
        Array containing the names of the equilibrium arrays.
    grid : numpy.ndarray(dtype=float, ndim=1)
        Array containing the base grid.
    grid_gauss : numpy.ndarray(dtype=float, ndim=1)
        Array containing the Gaussian grid.
    equilibria : dict
        Dictionary containing the equilibrium arrays.
    eigenvalues : numpy.ndarray(dtype=complex, ndim=1)
        Array containing the eigenvalues.
    continua : dict
        Dictionary containing the continuum regions.
    """
    def __init__(self, datfile):
        self.datfile = str(datfile)
        self._check_file_validity()
        with open(self.datfile, 'rb') as istream:
            self.header = get_header(istream)
        self._set_header_attributes()
        self._load_basic_data()
        self._calculate_continua()

    def _set_header_attributes(self):
        """
        Sets basic header data as instance attributes.
        """
        self.geometry = self.header['geometry']
        self.x_start = self.header['x_start']
        self.x_end = self.header['x_end']
        self.gridpts = self.header['gridpts']
        self.gauss_gridpts = self.header['gauss_gridpts']
        self.matrix_gridpts = self.header['matrix_gridpts']
        self.ef_gridpts = self.header['ef_gridpts']
        self.gamma = self.header['gamma']
        self.eq_type = self.header['eq_type']
        self.parameters = self.header['params']
        self.cgs = self.header['cgs']
        self.units = self.header['units']
        self.eq_names = self.header['equil_names']
        pylboLogger.info('gridpoints = {}'.format(self.gridpts))
        pylboLogger.info('geometry   = {}'.format(self.geometry))
        pylboLogger.info('grid start = {}'.format(self.x_start))
        pylboLogger.info('grid end   = {}'.format(self.x_end))

    def _load_basic_data(self):
        """
        Sets the base grid (:attr:`grid`), Gaussian grid (:attr:`grid_gauss`),
        equilibrium arrays (:attr:`equilibria`) and eigenvalues (:attr:`eigenvalues`).
        """
        with open(self.datfile, 'rb') as istream:
            self.grid = read_grid(istream, self.header)
            self.grid_gauss = read_grid_gauss(istream, self.header)
            self.equilibria = read_equilibrium_arrays(istream, self.header)
            self.eigenvalues = read_eigenvalues(istream, self.header)

    def _calculate_continua(self):
        """
        Calculates the various continua for this dataset and sets them
        as an instance attribute (:attr:`continua`), which contains the different regions
        as a dictionary.
        """
        wS_pos, wS_neg, wA_pos, wA_neg, wth = get_continuum_regions(self)
        self.continua = {'wS+': wS_pos, 'wS-': wS_neg, 'wA+': wA_pos,
                         'wA-': wA_neg, 'wth': wth}

    def get_sound_speed(self):
        """
        Calculates the sound speed based on the equilibrium arrays.

        Returns
        -------
        cs : numpy.ndarray(dtype=float, ndim=1)
            Sound speed at every grid point, defined as
            :math:`c_s = \\sqrt{\\frac{\\gamma p_0}{\\rho_0}}`.
        """
        pressure = self.equilibria['T0'] * self.equilibria['rho0']
        cs = np.sqrt(self.gamma * pressure / self.equilibria['rho0'])
        return cs

    def get_alfven_speed(self):
        """
        Calculates the Alfvén speed based on the equilibrium arrays.

        Returns
        -------
        vA : numpy.ndarray(dtype=float, ndim=1)
            Alfvén speed at every grid point, defined as
            :math:`v_A = \\sqrt{\\frac{B_0^2}{\\rho_0}}`.
        """
        B0 = np.sqrt(self.equilibria['B02']**2 + self.equilibria['B03']**2)
        vA = B0 / np.sqrt(self.equilibria['rho0'])
        return vA

    def get_k0_squared(self):
        """
        Calculates the squared wave number.

        Returns
        -------
        k0_sq : float
            The wave number squared, defined as
            :math:`k_0^2 = k_2^2 + k_3^2`.
        """
        k0_sq = self.parameters.get('k2')**2 + self.parameters.get('k3')**2
        return k0_sq

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
        if not isinstance(ev_guesses, np.ndarray):
            if isinstance(ev_guesses, (float, int, complex)):
                ev_guesses = [ev_guesses]
            ev_guesses = np.array(ev_guesses)
        idxs = np.empty(shape=len(ev_guesses), dtype=np.int)
        eigenvals = np.empty(shape=len(ev_guesses), dtype=np.complex)
        for i, guess in enumerate(ev_guesses):
            # get distance from guess to all eigenvalues
            distances = np.sqrt((self.eigenvalues.real - guess.real) ** 2
                                + (self.eigenvalues.imag - guess.imag) ** 2)
            # index of point with closest distance
            idx = distances.argmin()
            idxs[i] = idx
            eigenvals[i] = self.eigenvalues[idx]
        return idxs, eigenvals

    def get_matrix_B(self):
        """
        Retrieves the matrix B from the datfile.

        Returns
        -------
        matrix_B : numpy.ndarray(dtype=float, ndim=2)
            Matrix containing the B-matrix.

        Raises
        ------
        MatricesNotPresent
            If the matrices were not saved to the datfile.
        """
        if not self.header['matrices_written']:
            raise MatricesNotPresent(self.datfile)
        with open(self.datfile, 'rb') as istream:
            matrix_B = read_matrix_B(istream, self.header)
        return matrix_B

    def get_matrix_A(self):
        """
        Retrieves the matrix A from the datfile.

        Returns
        -------
        matrix_A : numpy.ndarray(dtype=complex, ndim=2)
            Matrix containing the A-matrix.

        Raises
        ------
        MatricesNotPresent
            If the matrices were not saved to the datfile.
        """
        if not self.header['matrices_written']:
            raise MatricesNotPresent(self.datfile)
        with open(self.datfile, 'rb') as istream:
            matrix_A = read_matrix_A(istream, self.header)
        return matrix_A

    def _check_file_validity(self):
        """
        Checks the file validity of a given datfile.

        Raises
        -------
        FileNotFoundError
            If the datfile can not be found.
        InvalidLegolasFile
            If the datfile is not a valid Legolas output file.
        """
        path_to_file = Path(self.datfile).resolve()
        if not path_to_file.is_file():
            raise FileNotFoundError(path_to_file)
        if not self.datfile.endswith('.dat'):
            raise InvalidLegolasFile(path_to_file)
