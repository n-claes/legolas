from utilities.parameters import DEFAULT_PARAMS
from utilities.tools import custom_enumerate, progressbar

import numpy as np
import f90nml
import subprocess
import os

if __name__ == '__main__':
    print('------------------------------------------')
    print('---- SCRIPT FOR ADVANCED LEGOLAS RUNS ----')
    print('------------------------------------------\n')

    equils = iter(DEFAULT_PARAMS.keys())
    for idx, eq in custom_enumerate(equils, start=1, step=3):
        print('{:03d} {:<30}{:03d} {:<30}{:03d} {}'.format(idx, eq, idx+1, next(equils, ""), idx+2, next(equils, "")))

    eq_nb = int(input("\nChoose equilibrium above: "))
    chosen_eq = list(DEFAULT_PARAMS.keys())[eq_nb-1]
    print('Selected equilibrium: {}'.format(chosen_eq))

    print("\nAvailable parameters:")
    params = iter(DEFAULT_PARAMS[chosen_eq].keys())
    for param in params:
        print('{:<30}{}'.format(param, next(params, "")))

    varied_param = None
    while varied_param not in list(DEFAULT_PARAMS[chosen_eq].keys()):
        varied_param = input('Which parameter should be varied?: ')
    lb = float(input('left bound: '))
    rb = float(input('right bound: '))

    nb_runs = int(input("\nHow many runs?: "))
    gridpts = int(input("Gridpoints for each run?: "))

    param_list = np.linspace(lb, rb, nb_runs)

    # generate parfiles and run legolas for each file. Do this from top directory.
    os.chdir('..')
    namelist = {'gridlist': {}, 'equilibriumlist': {}, 'savelist': {}, 'filelist': {}}
    print('\nRunning Legolas...')
    for run, param in enumerate(param_list, start=1):
        base_filename = 'multiruns/{:03d}_{}_{}_'.format(run, chosen_eq, varied_param)

        namelist['gridlist']['gridpoints'] = gridpts
        namelist['equilibriumlist']['equilibrium_type'] = chosen_eq
        namelist['savelist']['run_silent'] = True
        namelist['savelist']['show_results'] = False
        namelist['filelist']['savename_eigenvalues'] = base_filename + 'eigenvals'
        namelist['filelist']['savename_config'] = base_filename + 'config'

        # get default parameters for this equilibrium, add to namelist
        namelist['paramlist'] = DEFAULT_PARAMS[chosen_eq]

        # update parameter with changed value
        namelist['paramlist'].update({varied_param: param})

        # write the parfile, save to multiruns folder
        parfile_name = 'output/' + base_filename + 'parfile.par'
        f90nml.write(namelist, parfile_name, force=True)

        # run legolas
        progressbar(run, nb_runs, '{}/{}'.format(run, nb_runs))
        subprocess.check_call(['./legolas', '-i', parfile_name])
