import os
import subprocess
import signal
import multiprocessing
import f90nml
import psutil
import numpy as np
from pathlib import Path
from tqdm import tqdm
from .defaults import precoded_runs, \
    get_precoded_run, \
    LEGOLAS_DIR, \
    LEGOLAS_OUT
from .exceptions import DictNotEmpty, \
    InconsistentNumberOfRuns
from ..data_management.file_handler import select_files
from ..utilities.logger import pylboLogger

LEGOLAS_PAR = (LEGOLAS_OUT / 'parfiles').resolve()


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
    legolas_exec = (LEGOLAS_DIR / 'legolas').resolve()
    if not Path.is_file(legolas_exec):
        raise FileNotFoundError('Legolas executable not found in {}'.format(legolas_exec))

def custom_enumerate(iterable, start=0, step=1):
    """
    Does a custom enumeration with a given stepsize.

    Parameters
    ----------
    iterable : ~typing.Iterable
        The iterable to iterate over.
    start : int
        The starting value for enumerate.
    step : int
        The stepsize between enumerate values.

    Yields
    ------
    start : :class:`int`
        The current index in `iterable`, incremented with `step`.
    itr : :class:`~typing.Iterable`
        The corresponding entry of `iterable`.
    """
    for itr in iterable:
        yield start, itr
        start += step


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
    output_dir : str or Pathlike object
        The path to the output directory.

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
    def update_namelist(_key, _items):
        namelist.update({_key: {}})
        for _item in _items:
            _value = parfile_dict.pop(_item, None)
            if _value is not None:
                namelist[_key].update({_item: _value})
        if namelist[_key] == {}:
            namelist.pop(_key)

    if parfile_dict is None:
        print('No dictionary supplied. Available precoded runs:')
        av_runs = iter(precoded_runs.keys())
        for idx, pr in custom_enumerate(av_runs, start=1, step=3):
            print('{:03d} {:<30}{:03d} {:<30}{:03d} {}'.format(idx, pr, idx + 1, next(av_runs, "<empty>"),
                                                               idx + 2, next(av_runs, '<empty>')))
        pr_in = int(input('\nChoose precoded run: '))
        chosen_pr = list(precoded_runs.keys())[pr_in - 1]
        print('Selected run: {}'.format(chosen_pr))
        parfile_dict = get_precoded_run(chosen_pr)
    if output_dir is None:
        output_dir = LEGOLAS_PAR
    if isinstance(output_dir, str):
        output_dir = Path(output_dir).resolve()
    _check_directories(output_dir)
    nb_runs = parfile_dict.pop('number_of_runs', 1)
    namelist = {}
    update_namelist('gridlist', ['geometry', 'x_start', 'x_end', 'gridpoints',
                                 'mesh_accumulation', 'ev_1', 'ev_2', 'sigma_1',
                                 'sigma_2', 'force_r0'])
    update_namelist('equilibriumlist', ['equilibrium_type', 'boundary_type', 'use_defaults',
                                        'remove_spurious_eigenvalues', 'nb_spurious_eigenvalues'])
    update_namelist('savelist', ['write_matrices', 'write_eigenfunctions', 'show_results',
                                 'basename_datfile', 'basename_logfile', 'output_folder', 'logging_level'])
    update_namelist('physicslist', ['mhd_gamma', 'flow', 'radiative_cooling', 'ncool', 'cooling_curve',
                                    'external_gravity', 'thermal_conduction', 'use_fixed_tc_para', 'fixed_tc_para_value',
                                    'use_fixed_tc_perp', 'fixed_tc_perp_value', 'resistivity', 'use_fixed_resistivity',
                                    'fixed_eta_value'])
    update_namelist('unitslist', ['cgs_units', 'unit_density', 'unit_temperature', 'unit_magneticfield', 'unit_length'])
    if parfile_dict.get('parameters') is not None:
        namelist.update({'paramlist': parfile_dict.pop('parameters')})
    # we should have popped everything from the dictionary so it should be empty.
    # if it's not, something is wrong
    if not parfile_dict == {}:
        raise DictNotEmpty(parfile_dict)
    # explicitly turn off defaults if paramlist is supplied
    if namelist.get('paramlist') is not None:
        namelist['equilibriumlist'].update({'use_defaults': False})
    # create specific dictionary so multiple keys can be varied at the same time
    for major_key in namelist.keys():
        for key, item in namelist[major_key].items():
            if not isinstance(item, np.ndarray):
                namelist[major_key].update({key: [item] * nb_runs})
            else:
                if len(item) != nb_runs:
                    raise InconsistentNumberOfRuns(nb_runs, key, namelist[major_key])
    # generate parfiles
    parfiles = []
    for run in range(nb_runs):
        run_prepended = "{:04d}".format(run)
        if nb_runs == 1:
            run_prepended = ''
        parfile_dict = {key: {} for key in namelist.keys()}
        # create dictionary for single run
        for name, nl_dict in namelist.items():
            for key, item in nl_dict.items():
                parfile_dict[name].update({key: item[run]})
        if basename_parfile is None:
            basename_parfile = parfile_dict['equilibriumlist']['equilibrium_type']
        parfile_name = "{}{}.par".format(run_prepended, basename_parfile)
        basename_datfile = parfile_dict.get('savelist', {}).get('basename_datfile')
        if basename_datfile is None:
            basename_datfile = parfile_dict['equilibriumlist']['equilibrium_type']
        datfile_name = "{}{}".format(run_prepended, basename_datfile)
        parfile_dict['savelist'].update({'basename_datfile': datfile_name})
        basename_logfile = parfile_dict.get('savelist', {}).get('basename_logfile')
        if basename_logfile is not None:
            logfile_name = "{}{}".format(run_prepended, basename_logfile)
            parfile_dict['savelist'].update({'basename_logfile': logfile_name})
        # set paths, write parfile
        parfile_path = (output_dir / parfile_name).resolve()
        parfiles.append(str(parfile_path))
        f90nml.write(parfile_dict, parfile_path, force=True)
    pylboLogger.info('parfiles generated, saved to {}'.format(output_dir))
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
    cmd = ['./legolas', '-i', str(parfile)]
    return subprocess.call(cmd)


def _terminate_workers():
    """
    Terminates the multiprocessing workers after a forced interruption.
    Simply giving an interrupt terminates only the Python processes, but still
    keeps the Legolas calls running since those are subprocesses. This method first
    terminates all child processes (legolas), then the parents (workers).
    """
    pylboLogger.error('interrupting processes...')
    for process in multiprocessing.active_children():
        pid = process.pid
        pylboLogger.error('terminating PID: {} -- {}'.format(pid, process.name))
        parent = psutil.Process(pid)
        children = parent.children(recursive=True)
        for child in children:
            pylboLogger.error('terminating child process {} -- {}'.format(child.pid, child.name()))
            child.kill()
        gone, alive = psutil.wait_procs(children, timeout=2)
        for killed_proc in gone:
            pylboLogger.error('{}'.format(str(killed_proc)))
        parent.kill()
        parent.wait(timeout=2)
    pylboLogger.critical('all Legolas processes terminated.')


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
        pbar = tqdm(total=len(parfiles), unit='')
        pbar.set_description('Running Legolas [CPUs={}]'.format(nb_cpus))
        pool = multiprocessing.Pool(processes=nb_cpus, initializer=_init_worker)
        try:
            # start multiprocessing pool
            for parfile in parfiles:
                pool.apply_async(_activate_worker, args=(parfile,), callback=update_pbar)
            pool.close()
            pool.join()
            pbar.close()
        except KeyboardInterrupt:
            pbar.set_description('INTERRUPTED')
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
        pylboLogger.info('Parfiles removed.')
        # if parfile directory is empty, also remove that one
        try:
            Path.rmdir(LEGOLAS_PAR)
            pylboLogger.info('{} was empty, so also removed the directory.'.format(LEGOLAS_PAR))
        except OSError:
            pass
