import os
import subprocess
import signal
import multiprocessing
import f90nml
import psutil
from pathlib import Path
from tqdm import tqdm
from .toolbox import transform_to_list
from .defaults import (
    namelist_items,
    LEGOLAS_DIR,
    LEGOLAS_OUT,
)
from .exceptions import (
    DictNotEmpty,
    InconsistentNumberOfRuns,
)
from ..data_management.file_handler import select_files
from ..utilities.logger import pylboLogger

LEGOLAS_PAR = (LEGOLAS_OUT / "parfiles").resolve()


def _check_directories(output_dir):
    """
    Does checks on the existence of specified directories.
    If `output_dir` is not found, it is created.

    Parameters
    ----------
    output_dir : ~os.PathLike
        The output directory to save files.

    Raises
    ------
    NotADirectoryError
        If the Legolas directory could not be found.
    """
    if not Path.is_dir(LEGOLAS_DIR):
        raise NotADirectoryError(LEGOLAS_DIR)
    if not Path.is_dir(LEGOLAS_OUT):
        Path.mkdir(LEGOLAS_OUT)
    if not Path.is_dir(output_dir):
        Path.mkdir(output_dir)


def _check_executable():
    """
    Checks the Legolas executable.

    Raises
    -------
    FileNotFoundError
        If the Legolas executable could not be found.
    """
    legolas_exec = (LEGOLAS_DIR / "legolas").resolve()
    if not Path.is_file(legolas_exec):
        raise FileNotFoundError(f"Legolas executable not found in {legolas_exec}")


def generate_parfiles(parfile_dict=None, basename_parfile=None, output_dir=None):
    """
    Generates (series of) parfiles based on a given configuration dictionary. The files are saved
    using `basename_parfile` as a base name, and are written to the directory `output_dir`.
    If `output_dir` does not exist yet, it is created.
    Simply provide a dictionary with as keys the namelist items you want in the parfile, and
    this routine takes care of the rest.

    Parameters
    ----------
    parfile_dict : dict
        Dictionary containing the keys to be placed in the parfile.
    basename_parfile : str
        Base name for the parfile.
    output_dir : ~os.PathLike
        The path to the output directory.

    Notes
    -----
    This routine is quite flexible, meaning that you can specify both single-values and lists
    as dictionary entries. When there are lists present multiple parfiles will be generated.
    The only requirement is that list sizes are consistent across dictionary keys.

    Examples
    --------
    The example below will generate three parfiles, all three with an adiabatic homogeneous
    equilibrium state. The first one has 50 gridpoints and a wall at position 1, the second
    one has 100 gridpoints with a wall at position 5, and so on.

    >>> import pylbo
    >>> config = dict(equilibrium_type='adiabatic_homo',
    >>>               gridpoints=[50, 100, 150],
    >>>               x_end=[1, 5, 10])
    >>> output = 'path_to_my_folder'
    >>> myname = 'my_name'
    >>> parfiles = pylbo.generate_parfiles(parfile_dict=config,
    >>>                                    output_dir=output,
    >>>                                    basename_parfile=myname)

    Returns
    -------
    parfiles : list
        A list with the paths to the generated parfiles.

    Raises
    ------
    DictNotEmpty
        If there are still keys remaining in the given dictionary. If parfile generation is
        successful, we pop everything from a copy of the given dict.
    InconsistentNumberOfRuns
        When the `number_of_runs` key is inconsistent with array sizes in the parfile.
    """
    if parfile_dict is None:
        pylboLogger.error("no dictionary supplied, unable to generate parfile(s)!")
    if basename_parfile is None:
        basename_parfile = parfile_dict.get("equilibrium_type")
    if output_dir is None:
        output_dir = LEGOLAS_PAR
    output_dir = Path(output_dir).resolve()
    _check_directories(output_dir)
    nb_runs = parfile_dict.pop("number_of_runs", 1)

    # create namelist format
    namelist = {}
    for key, names in namelist_items.items():
        # create dictionary entry
        namelist.update({key: {}})
        # loop over items in this specific namelist item
        for name in names:
            obj = parfile_dict.pop(name, None)
            if obj is not None:
                obj = transform_to_list(obj)
                namelist[key].update({name: obj})
        # if it's still an empty dictionary, remove it
        if namelist[key] == {}:
            namelist.pop(key)
    # add parameters if present
    param_dict = parfile_dict.pop("parameters", {})
    if len(param_dict) != 0:
        namelist["equilibriumlist"].update({"use_defaults": [False]})
        for name, param in param_dict.items():
            obj = transform_to_list(param)
            param_dict.update({name: obj})
    namelist.update({"paramlist": param_dict})
    # at this point the original dictionary should be empty, if it's not
    # there is a problem
    if not len(parfile_dict) == 0:
        raise DictNotEmpty(parfile_dict)

    # now update namelist for number of runs consistency,
    # at this point all items in the namelist are lists
    for listkey, name in namelist.items():
        for key, item in name.items():
            if len(item) == 1:
                itemlist = item * nb_runs
            elif len(item) == nb_runs:
                itemlist = item
            else:
                raise InconsistentNumberOfRuns(nb_runs, key, namelist[listkey])
            namelist[listkey].update({key: itemlist})

    # generate parfiles
    parfiles = []
    for run in range(nb_runs):
        number = "{:04d}".format(run)
        if nb_runs == 1:
            number = ""
        # dictionary for this run
        run_dict = {key: {} for key in namelist.keys()}
        for listkey, name in namelist.items():
            for key, item in name.items():
                run_dict[listkey].update({key: item[run]})
        # make sure savelist is present
        try:
            run_dict["savelist"]
        except KeyError:
            run_dict.update({"savelist": {}})
        # parfile name
        parfile_name = f"{number}{basename_parfile}.par"
        # datfile name
        basename_datfile = run_dict.get("savelist").get(
            "basename_datfile", basename_parfile
        )
        datfile_name = f"{number}{basename_datfile}"  # no .dat extension here
        run_dict["savelist"].update({"basename_datfile": datfile_name})
        # logfile name
        basename_logfile = run_dict.get("savelist").get("basename_logfile", None)
        if basename_logfile is not None:
            logfile_name = f"{number}{basename_logfile}"  # no .log extension here
            run_dict["savelist"].update({"basename_logfile": logfile_name})
        # set paths, write parfile
        parfile_path = (output_dir / parfile_name).resolve()
        parfiles.append(str(parfile_path))
        f90nml.write(run_dict, parfile_path, force=True)
    pylboLogger.info("parfiles generated, saved to {}".format(output_dir))
    return parfiles


def _init_worker():
    """
    Worker initialisation for the multiprocessing module.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def _activate_worker(parfile):
    """
    Worker activation for the multiprocessing module.
    Calls the Legolas executable as a subprocess, with the parfile
    as argument.

    Parameters
    ----------
    parfile : str or ~os.PathLike
        The path to the parfile.

    Returns
    -------
    call : :func:`subprocess.call`
        A call to a subprocess to run Legolas.
    """
    cmd = ["./legolas", "-i", str(parfile)]
    return subprocess.call(cmd)


def _terminate_workers():
    """
    Terminates the multiprocessing workers after a forced interruption.
    Simply giving an interrupt terminates only the Python processes, but still
    keeps the Legolas calls running since those are subprocesses. This method first
    terminates all child processes (legolas), then the parents (workers).
    """
    pylboLogger.error("interrupting processes...")
    for process in multiprocessing.active_children():
        pid = process.pid
        pylboLogger.error(f"terminating PID: {pid} -- {process.name}")
        parent = psutil.Process(pid)
        children = parent.children(recursive=True)
        for child in children:
            pylboLogger.error(f"terminating child process {child.pid} -- {child.name}")
            child.kill()
        gone, alive = psutil.wait_procs(children, timeout=2)
        for killed_proc in gone:
            pylboLogger.error(f"{str(killed_proc)}")
        parent.kill()
        parent.wait(timeout=2)
    pylboLogger.critical("all Legolas processes terminated.")
    exit(1)


def run_legolas(parfiles=None, remove_parfiles=False, nb_cpus=1):
    """
    Runs the legolas executable for a given list of parfiles. If more than
    one parfile is passed, the runs are done using the multiprocessing module.
    This can be further controlled by the `nb_cpus` kwarg. If multiprocessing
    is done, a progressbar is printed. Every CPU will have a single Legolas
    instance associated with it.

    Parameters
    ----------
    parfiles : list or numpy.ndarray
        A list or array containing the paths to the parfiles.
    remove_parfiles : bool
        If `True`, removes the parfile after running Legolas. Will also remove
        the containing folder if it turns out to be empty after parfile removal.
    nb_cpus : int
        The amount of CPUs to use when running Legolas. If equal to 1, no
        multiprocessing is done.
    """

    def update_pbar(*args):
        pbar.update()

    _check_executable()
    if parfiles is None:
        parfiles = select_files()
    # original working directory
    owd = os.getcwd()
    # change to source directory
    os.chdir(LEGOLAS_DIR)
    # don't use multiprocessing if there is only one job
    if len(parfiles) == 1 or nb_cpus == 1:
        parfile = parfiles[0]
        pylboLogger.info("running Legolas...")
        try:
            _activate_worker(parfile)
        except KeyboardInterrupt:
            _terminate_workers()
            exit(1)
    else:
        # initialise progressbar and multiprocessing pool
        pbar = tqdm(total=len(parfiles), unit="")
        pbar.set_description(f"Running Legolas [CPUs={nb_cpus}]")
        pool = multiprocessing.Pool(processes=nb_cpus, initializer=_init_worker)
        try:
            # start multiprocessing pool
            for parfile in parfiles:
                pool.apply_async(
                    _activate_worker, args=(parfile,), callback=update_pbar
                )
            pool.close()
            pool.join()
            pbar.close()
        except KeyboardInterrupt:
            pbar.set_description("INTERRUPTED")
            pbar.update(len(parfiles))
            pbar.close()
            _terminate_workers()
            pool.terminate()
            pool.join()
            exit(1)
    pylboLogger.info("all runs completed")
    # change back to original directory
    os.chdir(owd)
    # remove parfiles if asked
    if remove_parfiles:
        for file in parfiles:
            os.remove(file)
        pylboLogger.info("Parfiles removed.")
        # if parfile directory is empty, also remove that one
        try:
            Path.rmdir(LEGOLAS_PAR)
            pylboLogger.info(f"{LEGOLAS_PAR} was empty, so also removed the directory.")
        except OSError:
            pass
