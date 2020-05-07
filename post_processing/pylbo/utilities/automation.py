import sys
import os
import subprocess
import f90nml
import numpy as np
from pathlib import Path
from .defaults import precoded_runs
from .exceptions import DictNotEmpty
from ..data_management.file_handler import select_files

SRC_DIR = Path(__file__).parents[3]
OUT_DIR = (SRC_DIR / 'output').resolve()
PAR_DIR = (OUT_DIR / 'parfiles').resolve()

def _check_directories():
    if not Path.is_dir(SRC_DIR):
        raise NotADirectoryError(SRC_DIR)
    if not Path.is_dir(OUT_DIR):
        Path.mkdir(OUT_DIR)
    if not Path.is_dir(PAR_DIR):
        Path.mkdir(PAR_DIR)

def _check_executable():
    legolas_exec = (SRC_DIR / 'legolas').resolve()
    if not Path.is_file(legolas_exec):
        raise FileNotFoundError('Legolas executable not found in {}'.format(legolas_exec))

def custom_enumerate(iterable, start=0, step=1):
    for itr in iterable:
        yield start, itr
        start += step

def progressbar(current, total, text=''):
    bar_length = 40
    progress = float(current) / float(total)
    if progress >= 1.:
        progress, status = 1, "\r\n"
    block = int(round(bar_length * progress))
    text = "\r[{}] {:.0f}% {}".format("#" * block + "-" * (bar_length - block), round(progress * 100, 0), text)
    sys.stdout.write(text)
    sys.stdout.flush()

def run_legolas(parfiles=None, remove_parfiles=False):
    if parfiles is None:
        parfiles = select_files()

    print('Running Legolas...')
    nb_runs = len(parfiles)
    orig_dir = os.getcwd()
    # change to source directory
    os.chdir(SRC_DIR)
    # run legolas
    for run in range(nb_runs):
        progressbar(run, nb_runs, '{}/{}'.format(run, nb_runs))
        subprocess.check_call(['./legolas', '-i', parfiles[run]])
    progressbar(nb_runs, nb_runs, '{}/{}'.format(nb_runs, nb_runs))
    # change back to original directory
    print('\nDone.')
    os.chdir(orig_dir)

    # remove parfiles if asked
    if remove_parfiles:
        for file in parfiles:
            os.remove(file)
        print('Parfiles removed.')
        # if parfile directory is empty, also remove that one
        try:
            Path.rmdir(PAR_DIR)
            print('{} was empty, so also removed the directory.'.format(PAR_DIR))
        except OSError:
            # don't remove if director is not empty
            pass

def generate_parfiles(filename, config_dict):
    def update_namelist(name, items):
        namelist.update({name: {}})
        for item in items:
            _value = config_dict.pop(item, None)
            if _value is not None:
                namelist[name].update({item: _value})

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
                                 'ev_1', 'ev_2', 'sigma_1', 'sigma_2'])
    # handle equilibriumlist
    update_namelist('equilibriumlist', ['equilibrium_type', 'boundary_type', 'use_defaults'])
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
        parfile_path = (PAR_DIR / parfile_name).resolve()
        parfiles.append(parfile_path)

        # write parfile
        f90nml.write(namelist, parfile_path, force=True)
    return parfiles

def generate_multirun(parfile_name=None, config_dict=None, remove_parfiles=False):
    print('No dictionary supplied. Available precoded runs:')
    if config_dict is None:
        av_runs = iter(precoded_runs.keys())
        for idx, pr in custom_enumerate(av_runs, start=1, step=3):
            print('{:03d} {:<30}{:03d} {:<30}{:03d} {}'.format(idx, pr, idx + 1, next(av_runs, "<empty>"),
                                                               idx + 2, next(av_runs, '<empty>')))
        pr_in = int(input('\nChoose precoded run: '))
        chosen_pr = list(precoded_runs.keys())[pr_in - 1]
        print('Selected run: {}'.format(chosen_pr))
        config_dict = precoded_runs.get(chosen_pr)
    parfiles = generate_parfiles(parfile_name, config_dict)
    run_legolas(parfiles, remove_parfiles)
