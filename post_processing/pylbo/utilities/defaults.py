import numpy as np
from pathlib import Path


#: Path to the Legolas directory.
LEGOLAS_DIR = Path(__file__).parents[3]
#: Path to the Legolas output directory.
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
        'basename_datfile': 'gravito_acoustic_NB10'
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
        'basename_datfile': 'gravito_mhd_beta1'
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
        'basename_datfile': 'gravito_mhd_beta50'
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
        'basename_datfile': 'interchange_modes_noshear'
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
        'basename_datfile': 'interchange_modes_03shear'
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
        'basename_datfile': 'constant_current_m-2'
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
        'basename_datfile': 'photospheric_flux_tube'
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
        'basename_datfile': 'coronal_flux_tube'
    }
}

namelist_items = {
    'gridlist': [
        'geometry',
        'x_start',
        'x_end',
        'gridpoints',
        'mesh_accumulation',
        'ev_1',
        'ev_2',
        'sigma_1',
        'sigma_2',
        'force_r0'
    ],
    'equilibriumlist': [
        'equilibrium_type',
        'boundary_type',
        'use_defaults',
        'remove_spurious_eigenvalues',
        'nb_spurious_eigenvalues'
    ],
    'savelist': [
        'write_matrices',
        'write_eigenfunctions',
        'show_results',
        'basename_datfile',
        'basename_logfile',
        'output_folder',
        'logging_level'
    ],
    'physicslist': [
        'mhd_gamma',
        'flow',
        'radiative_cooling',
        'ncool',
        'cooling_curve',
        'external_gravity',
        'thermal_conduction',
        'use_fixed_tc_para',
        'fixed_tc_para_value',
        'use_fixed_tc_perp',
        'fixed_tc_perp_value',
        'resistivity',
        'use_fixed_resistivity',
        'fixed_eta_value',
        'use_eta_dropoff',
        'dropoff_edge_dist',
        'dropoff_width'
    ],
    'unitslist': [
        'cgs_units',
        'unit_density',
        'unit_temperature',
        'unit_magneticfield',
        'unit_length'
    ]
}
