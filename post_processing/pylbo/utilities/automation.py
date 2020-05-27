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

LEGOLAS_PAR = (LEGOLAS_OUT / 'parfiles').resolve()

def _check_directories(output_dir):
    if not Path.is_dir(LEGOLAS_DIR):
        raise NotADirectoryError(LEGOLAS_DIR)
    if not Path.is_dir(LEGOLAS_OUT):
        Path.mkdir(LEGOLAS_OUT)
    if not Path.is_dir(output_dir):
        Path.mkdir(output_dir)

def _check_executable():
    legolas_exec = (LEGOLAS_DIR / 'legolas').resolve()
    if not Path.is_file(legolas_exec):
        raise FileNotFoundError('Legolas executable not found in {}'.format(legolas_exec))

def custom_enumerate(iterable, start=0, step=1):
    for itr in iterable:
        yield start, itr
        start += step

def generate_parfiles(parfile_dict=None, basename_parfile=None, output_dir=None):
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
        parfile_dict = {'gridlist': {}, 'equilibriumlist': {}, 'paramlist': {}, 'savelist': {}}
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
    print('Parfiles generated, saved to {}'.format(output_dir))
    return parfiles

def _init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def _activate_worker(parfile):
    cmd = ['./legolas', '-i', str(parfile)]
    return subprocess.call(cmd)

def run_legolas(parfiles=None, remove_parfiles=False, nb_cpus=1):
    def update_pbar(*args):
        pbar.update()

    _check_executable()
    if parfiles is None:
        parfiles = select_files()
    # original working directory
    owd = os.getcwd()
    # change to source directory
    os.chdir(LEGOLAS_DIR)
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
        # Note: simply calling pool.terminate() terminates ONLY the python processes,
        # but still keeps the legolas calls running since those are done using subprocess.
        # The following code first terminates all child processes (legolas), then the parents (workers)
        print('\nCaught KeyboardInterrupt:')
        for process in multiprocessing.active_children():
            pid = process.pid
            print('  Terminating ID: {} -- {}'.format(pid, process.name))
            parent = psutil.Process(pid)
            children = parent.children(recursive=True)
            for child in children:
                print('    Terminating child process {} -- {}'.format(child.pid, child.name()))
                child.kill()
            gone, alive = psutil.wait_procs(children, timeout=2)
            for killed_proc in gone:
                print('    {}'.format(str(killed_proc)))
            parent.kill()
            parent.wait(timeout=2)
        pool.terminate()
        print('All processes terminated.\n')
    # change back to original directory
    os.chdir(owd)
    # remove parfiles if asked
    if remove_parfiles:
        for file in parfiles:
            os.remove(file)
        print('Parfiles removed.')
        # if parfile directory is empty, also remove that one
        try:
            Path.rmdir(LEGOLAS_PAR)
            print('{} was empty, so also removed the directory.'.format(LEGOLAS_PAR))
        except OSError:
            pass
