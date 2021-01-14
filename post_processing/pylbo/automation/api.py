from pylbo.automation.generator import ParfileGenerator
from pylbo.automation.runner import LegolasRunner


def generate_parfiles(parfile_dict, basename=None, output_dir=None):
    """
    Generates parfiles based on a given configuration dictionary.
    No namelists have to be supplied, simply providing a dictionary with keys the
    namelist items you want in the parfile and this routine takes care of the rest.

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
    """
    pfgen = ParfileGenerator(parfile_dict, basename, output_dir)
    pfgen.create_namelist_from_dict()
    parfiles = pfgen.generate_parfiles()
    return parfiles


def run_legolas(parfiles, remove_parfiles=False, nb_cpus=1, executable=None):
    """
    Runs the legolas executable for a given list of parfiles. If more than one parfile
    is passed, the runs can be performed in parallel using the multiprocessing module.
    Parallelisation can be enabled by setting the `nb_cpus` kwarg to a number greater
    than one. Every CPU will have a single legolas executable subprocess associated
    with it.

    Parameters
    ----------
    parfiles : list or numpy.ndarray
        A list or array containing the names or paths to the parfiles.
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
    """
    runner = LegolasRunner(parfiles, remove_parfiles, nb_cpus, executable)
    runner.execute()
