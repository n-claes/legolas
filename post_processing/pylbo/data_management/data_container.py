import numpy as np
from pathlib import Path
from utilities.api import get_header, \
    read_grid, \
    read_grid_gauss, \
    read_eigenvalues, \
    read_equilibrium_arrays, \
    read_ef_grid, \
    read_eigenfunctions, \
    read_matrix_B, \
    read_matrix_A
from utilities.api import InvalidLegolasFile, \
    EigenfunctionsNotPresent, \
    MatricesNotPresent
from utilities.api import get_continuum_regions


class LegolasDataContainer:
    def __init__(self, datfile):
        self.datfile = str(datfile)
        self._check_file_validity()

        with open(self.datfile, 'rb') as istream:
            self.header = get_header(istream)
        self._set_header_attributes()
        self._load_basic_data()
        self._calculate_continua()

    def _set_header_attributes(self):
        # grid attributes
        self.geometry = self.header['geometry']
        self.x_start = self.header['x_start']
        self.x_end = self.header['x_end']
        self.gridpts = self.header['gridpts']
        self.gauss_gridpts = self.header['gauss_gridpts']
        self.matrix_gridpts = self.header['matrix_gridpts']
        self.ef_gridpts = self.header['ef_gridpts']
        # physics attributes
        self.gamma = self.header['gamma']
        self.eq_type = self.header['eq_type']
        self.parameters = self.header['params']
        self.cgs = self.header['cgs']
        self.units = self.header['units']
        self.eq_names = self.header['equil_names']

    def _load_basic_data(self):
        with open(self.datfile, 'rb') as istream:
            self.grid = read_grid(istream, self.header)
            self.grid_gauss = read_grid_gauss(istream, self.header)
            self.equilibria = read_equilibrium_arrays(istream, self.header)
            self.eigenvalues = read_eigenvalues(istream, self.header)

    def _calculate_continua(self):
        # calculate continuum regions
        wS_pos, wS_neg, wA_pos, wA_neg, wth = get_continuum_regions(self)
        self.continua = {'wS+': wS_pos, 'wS-': wS_neg, 'wA+': wA_pos,
                         'wA-': wA_neg, 'wth': wth}

    def get_sound_speed(self):
        pressure = self.equilibria['T0'] * self.equilibria['rho0']
        cs = np.sqrt(self.gamma * pressure / self.equilibria['rho0'])
        return cs

    def get_alfven_speed(self):
        B0 = np.sqrt(self.equilibria['B02']**2 + self.equilibria['B03']**2)
        vA = B0 / np.sqrt(self.equilibria['rho0'])
        return vA

    def get_k0_squared(self):
        k0_sq = self.parameters.get('k2')**2 + self.parameters.get('k3')**2
        return k0_sq

    def get_eigenfunctions(self):
        if not self.header['eigenfuncs_written']:
            raise EigenfunctionsNotPresent(self.datfile)
        with open(self.datfile, 'rb') as istream:
            ef_grid = read_ef_grid(istream, self.header)
            eigenfunctions = read_eigenfunctions(istream, self.header)
        return ef_grid, eigenfunctions

    def get_matrix_B(self):
        if not self.header['matrices_written']:
            raise MatricesNotPresent(self.datfile)
        with open(self.datfile, 'rb') as istream:
            matrix_B = read_matrix_B(istream, self.header)
        return matrix_B

    def get_matrix_A(self):
        if not self.header['matrices_written']:
            raise MatricesNotPresent(self.datfile)
        with open(self.datfile, 'rb') as istream:
            matrix_A = read_matrix_A(istream, self.header)
        return matrix_A

    def _check_file_validity(self):
        path_to_file = Path(self.datfile).resolve()
        if not path_to_file.is_file():
            raise FileNotFoundError(path_to_file)
        if not self.datfile.endswith('.dat'):
            raise InvalidLegolasFile(path_to_file)
