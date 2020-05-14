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
    select_files, \
    read_log_file
from .utilities.automation import \
    generate_parfiles, \
    run_legolas
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
from .utilities.defaults import \
    LEGOLAS_DIR, \
    LEGOLAS_OUT
from .visualisations.spectrum import \
    SingleSpectrum, \
    MultiSpectrum
from .visualisations.matrices import plot_matrices

def show():
    plt.show()
