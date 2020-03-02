from utilities.parameters import DEFAULT_PARAMS, PRECODED_MULTIRUNS
from utilities.tools import custom_enumerate, progressbar
from pathlib import Path

import numpy as np
import f90nml
import subprocess
import sys
import os

# define directories
FILE_DIR = Path(__file__).parent
SRC_DIR = (FILE_DIR / '..').resolve()
MULTI_DIR = (SRC_DIR / 'output/multiruns/').resolve()
if not Path(MULTI_DIR).is_dir():
    print('Output directory not found, make sure Legolas is compiled!')
    sys.exit(1)

# make subdirectories
Path.joinpath(MULTI_DIR, 'nml').mkdir(parents=True, exist_ok=True)
Path.joinpath(MULTI_DIR, 'par').mkdir(parents=True, exist_ok=True)
Path.joinpath(MULTI_DIR, 'dat').mkdir(parents=True, exist_ok=True)

NML_DIR = MULTI_DIR.name + '/nml/'
DAT_DIR = MULTI_DIR.name + '/dat/'


def _generate_parameter_dictionary(params, nb_runs):
    param_dict = {}

    # Update floats to an array containing the constant values.
    # This is done to support varying more than one parameter at once
    for key, value in zip(params.keys(), params.values()):
        if not isinstance(value, np.ndarray):
            param_dict.update({key: value * np.ones(nb_runs)})
        else:
            param_dict.update({key: value})

    return param_dict


def _generate_parfiles(gridpts, nb_runs, chosen_eq, params, precoded, **kwargs):
    # generate appropriate dictionary containing parameters to iterate over
    param_dict = _generate_parameter_dictionary(params, nb_runs)

    namelist = {'gridlist': {}, 'equilibriumlist': {}, 'savelist': {}, 'filelist': {}}
    if precoded is None:
        base_filename = chosen_eq
    else:
        base_filename = precoded

    print('\nRunning Legolas... ({} gridpts)'.format(gridpts))
    for run in range(1, nb_runs+1):
        filename = '{}_{:03d}_'.format(base_filename, run)

        namelist['gridlist']['gridpoints'] = gridpts
        namelist['equilibriumlist']['equilibrium_type'] = chosen_eq
        namelist['equilibriumlist']['use_defaults'] = kwargs.get('use_defaults', False)
        namelist['savelist']['run_silent'] = kwargs.get('run_silent', True)
        namelist['savelist']['show_results'] = kwargs.get('show_results', False)
        namelist['filelist']['savename_eigenvalues'] = DAT_DIR + filename + 'eigenvals'
        namelist['filelist']['savename_config'] = NML_DIR + filename + 'config'

        # get parameters for this iteration, add to namelist
        namelist_params = {}
        for key, value_array in zip(param_dict.keys(), param_dict.values()):
            namelist_params.update({key: value_array[run-1]})
        namelist['paramlist'] = namelist_params

        # write the parfile, save to multiruns folder
        parfile_name = (Path.joinpath(MULTI_DIR, 'par') / (filename + 'parfile.par')).resolve()
        f90nml.write(namelist, parfile_name, force=True)

        # run legolas
        progressbar(run, nb_runs, '{}/{}'.format(run, nb_runs))
        os.chdir(str(SRC_DIR))
        subprocess.check_call(['./legolas', '-i', parfile_name])
        os.chdir(str(FILE_DIR))



def _main():
    print('------------------------------------------')
    print('---- SCRIPT FOR ADVANCED LEGOLAS RUNS ----')
    print('------------------------------------------\n')

    use_precoded = input('Generate pre-coded multirun? ').lower() in ['yes', 'y']
    chosen_pr = None

    if use_precoded:
        precoded_runs = iter(PRECODED_MULTIRUNS.keys())
        for idx, pr in custom_enumerate(precoded_runs, start=1, step=3):
            print('{:03d} {:<30}{:03d} {:<30}{:03d} {}'.format(idx, pr, idx + 1, next(precoded_runs, ""),
                                                               idx + 2, next(precoded_runs, "")))
        pr_in = int(input("\nChoose pre-coded run above: "))
        chosen_pr = list(PRECODED_MULTIRUNS.keys())[pr_in - 1]
        print("Selected pre-coded run: {}".format(chosen_pr))

        current_eq = PRECODED_MULTIRUNS[chosen_pr]['equilibrium']
        # retrieve default parameters for this equilibrium
        params = DEFAULT_PARAMS[current_eq]
        # update parameters that have to be changed
        params.update(PRECODED_MULTIRUNS[chosen_pr]['parameters'])

        nb_runs = PRECODED_MULTIRUNS[chosen_pr]['nb_runs']
        gridpts = PRECODED_MULTIRUNS[chosen_pr]['gridpts']

    else:
        equils = iter(DEFAULT_PARAMS.keys())
        for idx, eq in custom_enumerate(equils, start=1, step=3):
            print('{:03d} {:<30}{:03d} {:<30}{:03d} {}'.format(idx, eq, idx + 1, next(equils, ""),
                                                               idx + 2, next(equils, "")))

        eq_nb = int(input("\nChoose equilibrium above: "))
        current_eq = list(DEFAULT_PARAMS.keys())[eq_nb - 1]
        print('Selected equilibrium: {}'.format(current_eq))

        # retrieve default parameters for this equilibrium
        params = DEFAULT_PARAMS[current_eq]

        print("\nAvailable parameters:")
        params_iter = iter(params.keys())
        for p in params_iter:
            print('{:<30}{}'.format(p, next(params_iter, "")))

        varied_param = None
        while varied_param not in list(params.keys()):
            varied_param = input('Which parameter should be varied?: ')
        lb = float(input('left bound: '))
        rb = float(input('right bound: '))

        nb_runs = int(input("\nHow many runs?: "))
        gridpts = int(input("Gridpoints for each run?: "))

        # update parameter in corresponding dict
        params.update({varied_param: np.linspace(lb, rb, nb_runs, endpoint=True)})

    _generate_parfiles(gridpts=gridpts, nb_runs=nb_runs, chosen_eq=current_eq, params=params, precoded=chosen_pr)


if __name__ == '__main__':
    _main()
