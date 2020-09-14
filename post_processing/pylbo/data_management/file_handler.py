import tkinter as tk
import numpy as np
from pathlib import Path
from tkinter import filedialog
from .data_container import LegolasDataContainer


def load(datfiles):
    """
    Loader for Legolas datfiles.

    Parameters
    ----------
    datfiles : str, ~os.PathLike, list of str, list of PathLike
        The path to the datfile, either as a :class:`str`
        or :class:`~os.PathLike` object,
        or as a list of :class:`str` or :class:`~os.PathLike` objects.

    Returns
    -------
    ds : (list of) :class:`~.pylbo.LegolasDataContainer`
        A (list of) :class:`~pylbo.LegolasDataContainer` instance(s),
        corresponding to the datfile(s) provided.
    """
    if isinstance(datfiles, list):
        datasets = []
        for datfile in datfiles:
            ds = LegolasDataContainer(datfile)
            datasets.append(ds)
        if len(datasets) == 1:
            datasets = datasets.pop()
        return datasets
    ds = LegolasDataContainer(datfiles)
    return ds


def read_log_file(file, sort=False):
    """
    Reads a Legolas log file.

    Parameters
    ----------
    file : str or ~os.PathLike
        The path to the logfile.
    sort : bool
        If `True`, sorts the eigenvalues in the logfile. Sorting is done first
        on the real part, then on the imaginary part.

    Returns
    -------
    eigenvalues : np.ndarray(dtype=complex, ndim=1)
        Array containing the eigenvalues from the logfile.

    Raises
    ------
    ValueError
        If the file is not a .log file.
    FileNotFoundError
        If the file could not be found.
    """
    if ".log" not in str(file):
        raise ValueError(f"This is not a .log file: {file}")
    filepath = Path(file).resolve()
    if not filepath.is_file():
        raise FileNotFoundError(filepath)
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
