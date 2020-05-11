import numpy as np
import copy
from pathlib import Path
from .exceptions import UnknownPrecodedRun

LEGOLAS_DIR = Path(__file__).parents[3]
LEGOLAS_OUT = (LEGOLAS_DIR / 'output').resolve()

precoded_runs = {
    'gravito_acoustic_NB10': {
        'geometry': 'Cartesian',
        'x_start': 0,
        'x_end': 1,
        'equilibrium_type': 'gravito_acoustic',
        'gridpoints': 31,
        'number_of_runs': 40,
        'parameters': {
            'k2': np.linspace(0, np.sqrt(250), 40),
            'k3': np.linspace(0, np.sqrt(250), 40),
            'cte_p0': 1,
            'g': 0.5,
            'alpha': 20.42
        },
        'savename_datfile': 'gravito_acoustic_NB10'
    },
    'gravito_mhd_beta1': {
        'geometry': 'Cartesian',
        'x_start': 0,
        'x_end': 1,
        'equilibrium_type': 'gravito_mhd',
        'gridpoints': 31,
        'number_of_runs': 40,
        'parameters': {
            'k2': np.linspace(0, np.sqrt(250), 40),
            'k3': np.linspace(0, np.sqrt(250), 40),
            'cte_p0': 0.5,
            'g': 20,
            'alpha': 20
        },
        'savename_datfile': 'gravito_mhd_beta1'
    },
    'gravito_mhd_beta50': {
        'geometry': 'Cartesian',
        'x_start': 0,
        'x_end': 1,
        'equilibrium_type': 'gravito_mhd',
        'gridpoints': 31,
        'number_of_runs': 40,
        'parameters': {
            'k2': np.linspace(0, np.sqrt(250), 40),
            'k3': np.linspace(0, np.sqrt(250), 40),
            'cte_p0': 25,
            'g': 20,
            'alpha': 20
        },
        'savename_datfile': 'gravito_mhd_beta50'
    },
    'interchange_modes_noshear': {
        'geometry': 'Cartesian',
        'x_start': 0,
        'x_end': 1,
        'equilibrium_type': 'interchange_modes',
        'gridpoints': 31,
        'number_of_runs': 40,
        'parameters': {
            'k2': np.pi * np.sin(np.linspace(0, np.pi, 40)),
            'k3': np.pi * np.cos(np.linspace(0, np.pi, 40)),
            'g': 0.5,
            'cte_p0': 0.25,
            'lambda': 0,
            'alpha': 20.0
        },
        'savename_datfile': 'interchange_modes_noshear'
    },
    'interchange_modes_0.3shear': {
        'geometry': 'Cartesian',
        'x_start': 0,
        'x_end': 1,
        'equilibrium_type': 'interchange_modes',
        'gridpoints': 31,
        'number_of_runs': 40,
        'parameters': {
            'k2': np.pi * np.sin(np.linspace(0, np.pi, 40)),
            'k3': np.pi * np.cos(np.linspace(0, np.pi, 40)),
            'g': 0.5,
            'cte_p0': 0.25,
            'lambda': 0.3,
            'alpha': 20.0
        },
        'savename_datfile': 'interchange_modes_03shear'
    },
    'constant_current_m-2': {
        'geometry': 'cylindrical',
        'x_start': 0,
        'x_end': 1,
        'equilibrium_type': 'constant_current_tokamak',
        'gridpoints': 51,
        'number_of_runs': 39,
        'parameters': {
            'k2': -2.0,
            'k3': 0.2,
            'j0': (2.0 * 0.2) / np.linspace(1.9, 2.1, 39)
        },
        'savename_datfile': 'constant_current_m-2'
    },
    'photospheric_flux_tube': {
        'geometry': 'cylindrical',
        'x_start': 0,
        'x_end': 2,
        'equilibrium_type': 'photospheric_flux_tube',
        'gridpoints': 51,
        'number_of_runs': 25,
        'parameters': {
            'k2': 0.0,
            'k3': np.linspace(0.1, 6.2, 25),
            'cte_rho0': 1.0,
            'cte_p0': 1.0,
            'r0': 1.0
        },
        'savename_datfile': 'photospheric_flux_tube'
    },
    'coronal_flux_tube': {
        'geometry': 'cylindrical',
        'x_start': 0,
        'x_end': 2,
        'equilibrium_type': 'coronal_flux_tube',
        'gridpoints': 51,
        'number_of_runs': 25,
        'parameters': {
            'k2': 0.0,
            'k3': np.linspace(0.1, 6.2, 25),
            'cte_rho0': 1.0,
            'cte_p0': 1.0,
            'r0': 1.0
        },
        'savename_datfile': 'coronal_flux_tube'
    }
}

def get_precoded_run(name, gridpoints=None, number_of_runs=None, savename_datfile=None):
    try:
        selected_run = copy.deepcopy(precoded_runs[name])
    except KeyError:
        raise UnknownPrecodedRun(name, precoded_runs.keys())
    for key, value in zip(['gridpoints', 'number_of_runs', 'savename_datfile'],
                          [gridpoints, number_of_runs, savename_datfile]):
        if value is not None:
            selected_run.update({key: value})
    return selected_run
