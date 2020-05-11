import os
import subprocess
import signal
import multiprocessing
import f90nml
import copy
import psutil
import numpy as np
from pathlib import Path
from tqdm import tqdm
from .defaults import precoded_runs, \
    get_precoded_run, \
    LEGOLAS_DIR, \
    LEGOLAS_OUT
from .exceptions import DictNotEmpty
from ..data_management.file_handler import select_files

LEGOLAS_PAR = (LEGOLAS_OUT / 'parfiles').resolve()

def _check_directories():
    if not Path.is_dir(LEGOLAS_DIR):
        raise NotADirectoryError(LEGOLAS_DIR)
    if not Path.is_dir(LEGOLAS_OUT):
        Path.mkdir(LEGOLAS_OUT)
    if not Path.is_dir(LEGOLAS_PAR):
        Path.mkdir(LEGOLAS_PAR)

def _check_executable():
    legolas_exec = (LEGOLAS_DIR / 'legolas').resolve()
    if not Path.is_file(legolas_exec):
        raise FileNotFoundError('Legolas executable not found in {}'.format(legolas_exec))

def custom_enumerate(iterable, start=0, step=1):
    for itr in iterable:
        yield start, itr
        start += step

def _init_work():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def _do_work(parfile):
    cmd = ['./legolas', '-i', str(parfile)]
    return subprocess.call(cmd)

def run_legolas(parfiles=None, remove_parfiles=False, cpus_to_use=None):
    def update_pbar(*args):
        pbar.update()

    if parfiles is None:
        parfiles = select_files()
    if cpus_to_use is None:
        cpus_to_use = int(multiprocessing.cpu_count() / 2)
    # save current working directory
    orig_dir = os.getcwd()
    # change to source directory
    os.chdir(LEGOLAS_DIR)
    # initialise progressbar and multiprocessing pool
    pbar = tqdm(total=len(parfiles), unit='')
    pbar.set_description('Running Legolas [CPUs={}]'.format(cpus_to_use))
    pool = multiprocessing.Pool(processes=cpus_to_use, initializer=_init_work)
    try:
        # start multiprocessing pool
        for parfile in parfiles:
            pool.apply_async(_do_work, args=(parfile,), callback=update_pbar)
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
    os.chdir(orig_dir)
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
            # don't remove if director is not empty
            pass

def generate_parfiles(config_dict, filename=None):
    def update_namelist(name, items):
        namelist.update({name: {}})
        for item in items:
            _value = config_dict.pop(item, None)
            if _value is not None:
                namelist[name].update({item: _value})

    config_dict = copy.deepcopy(config_dict)
    _check_directories()
    if filename is None:
        filename = config_dict.get('equilibrium_type')
    parfiles = []
    nb_runs = config_dict.pop('number_of_runs')

    # create specific dictionary for the parameters, so multiple parameters
    # can be varied for at the same time.
    params = config_dict.pop('parameters')
    for key, value in params.items():
        if not isinstance(value, np.ndarray):
            params.update({key: value * np.ones(nb_runs)})
        else:
            params.update({key: value})

    namelist = {}
    # handle gridlist
    update_namelist('gridlist', ['geometry', 'x_start', 'x_end',
                                 'gridpoints', 'mesh_accumulation',
                                 'ev_1', 'ev_2', 'sigma_1', 'sigma_2', 'force_r0'])
    # handle equilibriumlist
    update_namelist('equilibriumlist', ['equilibrium_type', 'boundary_type', 'use_defaults',
                                        'remove_spurious_eigenvalues', 'nb_spurious_eigenvalues'])
    # if not explicitly overridden don't use defaults
    if namelist['equilibriumlist'].get('use_defaults', None) is None:
        namelist['equilibriumlist'].update({'use_defaults': False})
    # handle savelist
    update_namelist('savelist', ['write_matrices', 'write_eigenfunctions',
                                 'show_results', 'savename_datfile', 'logging_level'])
    # if not explicitly overridden run silent
    if namelist['savelist'].get('logging_level') is None:
        namelist['savelist'].update({'logging_level': 0})
    # if not explicitly overriden do not show results
    if namelist['savelist'].get('show_results') is None:
        namelist['savelist'].update({'show_results': False})

    # we popped everything from the available dictionary, so it should be
    # empty by now. If not, something is wrong
    if not len(config_dict) == 0:
        raise DictNotEmpty(config_dict)

    datfile_basename = namelist['savelist'].get('savename_datfile', None)
    if datfile_basename is None:
        datfile_basename = filename
    # handle parameters, save parfiles
    namelist.update({'paramlist': {}})
    for run in range(nb_runs):
        for key, value_array in params.items():
            namelist['paramlist'].update({key: value_array[run]})
        parfile_name = '{}_{:04d}.par'.format(filename, run)
        datfile_name = '{}_{:04d}'.format(datfile_basename, run)
        namelist['savelist'].update({'savename_datfile': datfile_name})
        parfile_path = (LEGOLAS_PAR / parfile_name).resolve()
        parfiles.append(parfile_path)

        # write parfile
        f90nml.write(namelist, parfile_path, force=True)
    return parfiles

def generate_multirun(config_dict=None, parfile_name=None, remove_parfiles=False, gridpoints=None, datfile_name=None,
                      cpus_to_use=None):
    if config_dict is None:
        print('No dictionary supplied. Available precoded runs:')
        av_runs = iter(precoded_runs.keys())
        for idx, pr in custom_enumerate(av_runs, start=1, step=3):
            print('{:03d} {:<30}{:03d} {:<30}{:03d} {}'.format(idx, pr, idx + 1, next(av_runs, "<empty>"),
                                                               idx + 2, next(av_runs, '<empty>')))
        pr_in = int(input('\nChoose precoded run: '))
        chosen_pr = list(precoded_runs.keys())[pr_in - 1]
        print('Selected run: {}'.format(chosen_pr))
        config_dict = get_precoded_run(chosen_pr, gridpoints=gridpoints, savename_datfile=datfile_name)
    print('Running with {} gridpoints'.format(config_dict['gridpoints']))
    parfiles = generate_parfiles(config_dict=config_dict, filename=parfile_name)
    run_legolas(parfiles, remove_parfiles, cpus_to_use=cpus_to_use)
