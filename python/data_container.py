from utilities.read_data import read_stream_data
from pathlib import Path

import numpy as np

class LEGOLASDataContainer(list):
    def __init__(self, namelist_array):
        super().__init__()

        for namelist in namelist_array:
            self.append(_SingleDataContainer(namelist))


class _SingleDataContainer:
    def __init__(self, namelist):
        # source directory is two directories up
        self.src_dir = Path(__file__).parent.parent
        self.namelist = namelist

        self.gridpts = self.namelist['gridlist'].get('gridpoints')
        self.mat_gridpts = self.namelist['gridlist'].get('matrix_gridpts')
        self.ef_gridpts = self.namelist['gridlist'].get('ef_gridpts')
        self.current_eq = self.namelist['equilibriumlist'].get('current_equilibrium')

        self.output_folder = self.namelist['filelist'].get('output_folder')
        self.file_ext = self.namelist['filelist'].get('file_extension')

        self.fname_grid = self.output_folder + self.namelist['filelist'].get('savename_efgrid') + self.file_ext
        self.fname_matA = self.output_folder + self.namelist['filelist'].get('savename_matrix') + '_A' + self.file_ext
        self.fname_matB = self.output_folder + self.namelist['filelist'].get('savename_matrix') + '_B' + self.file_ext
        self.fname_w = self.output_folder + self.namelist['filelist'].get('savename_eigenvalues') + self.file_ext

        self.show_mats = self.namelist['savelist'].get('show_matrices')
        self.show_eigenfuncs = self.namelist['savelist'].get('show_eigenfunctions')
        self.ef_list = np.asarray(['rho', 'v1', 'v2', 'v3', 'T', 'a1', 'a2', 'a3'])

        self.params = self.namelist['paramlist']

        self.omegas = None
        self.matA = None
        self.matB = None
        self.grid = None
        self.eigenfuncs = None

        self.read_omegas()
        if self.show_mats:
            self.read_matrices()
        if self.show_eigenfuncs:
            self.read_eigenfuncs()


    def read_omegas(self):
        fname = (self.src_dir / self.fname_w).resolve()
        self.omegas = read_stream_data(fname, content_type=np.complex, rows=1, cols=self.mat_gridpts)[0]


    def read_matrices(self):
        fname = (self.src_dir / self.fname_matA).resolve()
        self.matA = read_stream_data(fname, content_type=np.complex, rows=self.mat_gridpts, cols=self.mat_gridpts)
        fname = (self.src_dir / self.fname_matB).resolve()
        self.matB = read_stream_data(fname, content_type=np.float64, rows=self.mat_gridpts, cols=self.mat_gridpts)


    def read_eigenfuncs(self):
        fname = (self.src_dir / self.fname_grid).resolve()
        self.grid = read_stream_data(fname, content_type=np.float64, rows=1, cols=self.ef_gridpts)[0]

        eigenfuncs = {}
        for i in range(len(self.ef_list)):
            varname = self.ef_list[i]
            ef_name = self.output_folder + 'eigenfunctions/{}_{}_eigenfunction'.format(i + 1, varname) + self.file_ext
            fname = (self.src_dir / ef_name).resolve()
            eigenfuncs[varname] = read_stream_data(fname, content_type=np.complex,
                                                   rows=self.mat_gridpts, cols=self.ef_gridpts)
        self.eigenfuncs = eigenfuncs
