from pylbo.automation.generator import ParfileGenerator
from pylbo.automation.runner import LegolasRunner


def generate_parfiles(parfile_dict, basename=None, output_dir=None, subdir=True):
    """
    Generates parfiles based on a given configuration dictionary.
    The separate namelists do not have to be taken into account, and a normal
    dictionary should be supplied where the keys correspond to the namelist
    items that are required. Typechecking is done automatically during parfile
    generation.

    Parameters
    ----------
    parfile_dict : dict
        Dictionary containing the keys to be placed in the parfile.
    basename : str
        The basename for the parfile, the `.par` suffix is added automatically and is
        not needed. If multiple parfiles are generated, these
        will be prepended by a 4-digit number (e.g. 0003myparfile.par).
        If not provided, the basename will default to `parfile`.
    output_dir : str, ~os.PathLike
        Output directory where the parfiles are saved, defaults to the current
        working directory if not specified. A subdirectory called `parfiles` will be
        created in which the parfiles will be saved.
    subdir : boolean
        If `.true.` (default), creates a subdirectory `parfiles` in the output folder.

    Notes
    -----
    This routine is quite flexible and specifically designed for parametric studies.
    You can specify both single values and list-like items as dictionary items.
    This routine will automatically generate multiple parfiles if lists/numpy arrays
    are present.

    Returns
    -------
    parfiles : list
        A list with the paths to the parfiles that were generated. This list can be
        passed immediately to :func:`pylbo.run_legolas`.

    Examples
    --------
    This will generate a single parfile in a subdirectory `parfile` of the
    current working directory.

    >>> import pylbo
    >>> config = {
    >>>    "geometry": "Cartesian",
    >>>    "x_start": 0,
    >>>    "x_end": 1,
    >>>    "equilibrium_type": "resistive_homo",
    >>>    "gridpoints": 100,
    >>>    "write_eigenfunctions": True,
    >>>    "basename_datfile": "my_run",
    >>>    "output_folder": "output",
    >>> }
    >>> parfile = pylbo.generate_parfiles(config)

    This will generate 10 parfiles in the directory `my_parfiles` relative to
    the current working directory. The first parfile will have `x_end = 1.0` and 100
    gridpoints, the second one will have `x_end = 2.0` and 150 gridpoints, etc.

    >>> import pylbo
    >>> import numpy as np
    >>> config = {
    >>>    "geometry": "Cartesian",
    >>>    "x_start": 0,
    >>>    "x_end": np.arange(1, 11)
    >>>    "number_of_runs": 10
    >>>    "equilibrium_type": "resistive_homo",
    >>>    "gridpoints": np.arange(100, 600, 50),
    >>>    "write_eigenfunctions": True,
    >>>    "basename_datfile": "my_run",
    >>>    "output_folder": "output",
    >>> }
    >>> parfile_list = pylbo.generate_parfiles(config, output_dir="my_parfiles")
    """
    pfgen = ParfileGenerator(parfile_dict, basename, output_dir, subdir)
    pfgen.create_namelist_from_dict()
    parfiles = pfgen.generate_parfiles()
    return parfiles


def run_legolas(parfiles, remove_parfiles=False, nb_cpus=1, executable=None):
    """
    Runs the legolas executable for a given list of parfiles. If more than one parfile
    is passed, the runs can be performed in parallel using the multiprocessing module.
    Parallelisation can be enabled by setting the `nb_cpus` kwarg to a number greater
    than one, and is disabled by default.
    Every CPU will have a single legolas executable subprocess associated
    with it.

    Parameters
    ----------
    parfiles : str, list, numpy.ndarray
        A string, list or array containing the paths to the parfile(s).
        Accepts the output of :func:`pylbo.generate_parfiles`.
    remove_parfiles : bool
        If `True`, removes the parfiles after running Legolas. This will also remove
        the containing folder if it turns out to be empty after the parfiles are
        removed. If there are other files still in the folder it remains untouched.
    nb_cpus : int
        The number of CPUs to use when running Legolas. If equal to 1 then
        parallelisation is disabled. Defaults to the maximum number of CPUs available
        if a number larger than the available number is specified.
    executable : str, ~os.PathLike
        The path to the legolas executable. If not specified, defaults to the
        standard one in the legolas home directory.

    Notes
    -----
    If multiprocessing is enabled, it is usually a good idea to have the number
    of runs requested divisible by the number of CPUs that are available. For example,
    if 24 runs are requested it is good practice to use either 2, 4, 6 or 8 CPUs,
    and avoid numbers like 3, 5 and 7.

    Examples
    --------
    The example below will run a list of parfiles using using a local legolas
    executable from the current directory, on 4 CPU's.

    >>> import pylbo
    >>> from pathlib import Path
    >>> files = sorted(Path("parfiles").glob("*.dat"))
    >>> pylbo.run_legolas(files, nb_cpus=4, executable="legolas")
    """
    runner = LegolasRunner(parfiles, remove_parfiles, nb_cpus, executable)
    runner.execute()
