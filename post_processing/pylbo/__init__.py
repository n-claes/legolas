import matplotlib.pyplot as plt
from . import \
    data_management, \
    testing, \
    utilities, \
    visualisations
from .data_management.data_container import \
    LegolasDataContainer
from .data_management.file_handler import \
    load, \
    select_files
from .utilities.automation import \
    run_legolas, \
    generate_parfiles, \
    generate_multirun
from .utilities.datfile_utils import \
    get_header, \
    read_grid, \
    read_grid_gauss, \
    read_eigenvalues, \
    read_equilibrium_arrays, \
    read_ef_grid, \
    read_eigenfunctions, \
    read_matrix_B, \
    read_matrix_A
from .visualisations.spectrum import \
    SingleSpectrum, \
    MultiSpectrum
from .visualisations.matrices import plot_matrices

from pathlib import Path

LEGOLAS_SRC = Path(__file__).parents[2].resolve()
LEGOLAS_OUT = Path(LEGOLAS_SRC / 'output').resolve()

def show():
    plt.show()
