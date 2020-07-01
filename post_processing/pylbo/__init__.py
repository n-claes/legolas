import matplotlib.pyplot as plt
from ._version import __version__
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
    read_eigenvalues
from .utilities.defaults import \
    LEGOLAS_DIR, \
    LEGOLAS_OUT
from .visualisations.spectrum import \
    SingleSpectrum, \
    MultiSpectrum
from .visualisations.matrices import \
    plot_matrices
from .visualisations.eigenfunctions import \
    EigenfunctionHandler
from .utilities.logger import \
    set_loglevel, \
    disable_logging

utilities.logger.init_logger()

def show():
    plt.show()
