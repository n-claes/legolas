import numpy as np
import tkinter as tk
from pathlib import Path
from tkinter import filedialog
from pylbo.exceptions import InvalidLegolasFile
from pylbo.data_containers import LegolasDataSet, LegolasDataSeries
from pylbo.utilities.toolbox import transform_to_list
from pylbo.utilities.logger import pylboLogger


def _validate_file(file):
    """
    Checks the file validity of a given logfile or datfile.

    Parameters
    ----------
    file : str, ~os.PathLike
        The path to the datfile, either as a :class:`str` or
        :class:`~os.PathLike` object.

    Raises
    -------
    FileNotFoundError
        If the datfile can not be found.
    InvalidLegolasFile
        If the datfile is not a valid Legolas output file.
    """
    path_to_file = Path(file).resolve()
    if not path_to_file.is_file():
        raise FileNotFoundError(path_to_file)
    if path_to_file.suffix not in (".dat", ".log"):
        raise InvalidLegolasFile(path_to_file)


def load(datfile):
    """
    Loads a single Legolas datfile.

    Parameters
    ----------
    datfile : str, ~os.PathLike
        Path to the datfile.

    Returns
    -------
    ds : :class:`~.pylbo.data_containers.DataSet`
        A dataset instance for the given datfile.
    """
    _validate_file(datfile)
    ds = LegolasDataSet(datfile)
    pylboLogger.info(f"Legolas version  : {ds.legolas_version}")
    pylboLogger.info(f"file loaded      : {ds.datfile.parent} -- {ds.datfile.name}")
    pylboLogger.info(f"gridpoints       : {ds.gridpts}")
    pylboLogger.info(f"geometry         : {ds.geometry} in {ds.x_start, ds.x_end}")
    pylboLogger.info(f"equilibrium      : {ds.eq_type}")
    if ds.header["matrices_written"]:
        pylboLogger.info("matrices present in datfile")
    if ds.header["eigenfuncs_written"]:
        pylboLogger.info("eigenfunctions present in datfile")
    print("-"*100)
    return ds


def load_series(datfiles):
    """
    Loads multiple Legolas datfiles.

    Parameters
    ----------
    datfiles : list of str, list of ~os.PathLike
        Paths to the datfiles that should be loaded, in list/array form.

    Returns
    -------
    series : :class:`~pylbo.data_containers.DataSeries`
        A dataseries instance for the given datfiles.
    """
    transform_to_list(datfiles)
    for datfile in datfiles:
        _validate_file(datfile)
    series = LegolasDataSeries(datfiles)
    return series


def load_logfile(logfile, sort=False):
    """
    Reads a Legolas log file.

    Parameters
    ----------
    logfile : str or ~os.PathLike
       The path to the logfile.
    sort : bool
       If `True`, sorts the eigenvalues in the logfile. Sorting is done first
       on the real part, then on the imaginary part.

    Returns
    -------
    eigenvalues : np.ndarray(dtype=complex, ndim=1)
       Array containing the eigenvalues from the logfile.
    """
    _validate_file(logfile)
    filepath = Path(logfile).resolve()
    eigenvalues = []
    with open(filepath, "r") as logfile:
        for line in logfile:
            line = line.strip().split(",")
            x, y = line
            eigenvalues.append(complex(float(x), float(y)))
    eigenvalues = np.asarray(eigenvalues)
    if sort:
        eigenvalues = np.sort(eigenvalues)
    return eigenvalues


def select_files():
    """
    Opens an interactive window to select files to open.
    Requires a visual interface.

    Returns
    -------
    files : list
        A list containing the paths to the files selected.
    """
    root = tk.Tk()
    root.withdraw()
    root.lift()
    root.focus_set()
    files = list(
        filedialog.askopenfilenames(parent=root, title="Select Legolas file(s)")
    )
    root.destroy()
    if not files:
        exit()
    return files
