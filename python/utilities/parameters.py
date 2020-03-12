import numpy as np

DEFAULT_PARAMS = {
    'adiabatic_homo': {
        'k2': np.pi,
        'k3': np.pi
    },
    'constant_current_tokamak': {
        'k2': np.pi,
        'k3': np.pi,
        'j0': 0.125
    },
    'discrete_alfven': {
        'k2': 1.0,
        'k3': 0.05,
        'nu': 2.0,
        'j0': 0.125,
        'r0': 0.2
    },
    'flow_driven_instabilities': {
        'k2': 0.0,
        'k3': 0.0,
        'g': 100.0,
        'delta': -5.0,
        'alpha': -np.pi,
        'theta': 0.0,
        'tau': 4.0,
        'p1': 0.5*np.pi,
        'p2': 1.0,
        'p3': 2.0,
        'p4': 1.0
    },
    'gravito_acoustic': {
        'k2': np.pi,
        'k3': np.pi,
        'cte_p0': 1.0,
        'g': 0.5,
        'alpha': 20.42
    },
    'gravito_mhd': {
        'k2': np.pi,
        'k3': np.pi,
        'cte_p0': 0.5,
        'g': 0.5,
        'alpha': 20.0
    },
    'ideal_quasimodes': {
        'k2': 2.0,
        'j0': 0.5,
        'r0': 10*np.pi,
        'p1': 1.0
    },
    'interchange_modes': {
        'k2': np.pi,
        'k3': np.pi,
        'g': 0.5,
        'cte_p0': 1.0,
        'lambda': 0.3,
        'alpha': 20.0
    },
    'interface_modes': {
        'k2': 0.0,
        'k3': 1.0,
        'cte_rho0': 1.0,
        'cte_T0': 1.0,
        'p1': 10.0,
        'p2': 10.0,
        'p3': 10.0
    },
    'internal_kink': {
        'k2': 1.0,
        'cte_rho0': 1.0,
        'cte_v03': 1.0,
        'cte_p0': 3.0
    },
    'kh': {
        'k2': 10.0,
        'k3': 0.0,
        'cte_p0': 3.6,
        'cte_v02': 1.67,
        'cte_v03': 0.0,
        'p1': 0.05
    },
    'kh_cd': {
        'k2': -1.0,
        'rc': 0.5,
        'V': 1.63,
        'cte_p0': 1.0,
        'Bz0': 0.25
    },
    'magneto_rotational': {
        'k2': 1.0,
        'k3': 70.0,
        'Bth0': 0.01,
        'Bz0': 0.01,
        'p1': 0.01,
        'p2': 1.0
    },
    'nonuniform_conduction': {
        'k2': 0.0,
        'k3': 1.0,
        'beta': 0.25
    },
    'resistive_homo': {
        'k2': 1.0,
        'k3': 1.0,
        'beta': 0.25
    },
    'resistive_tearing': {
        'k2': 0.49,
        'k3': 0.0,
        'alpha': 4.73884,
        'beta': 0.15
    },
    'resistive_tearing_flow': {
        'k2': 1.5,
        'k3': 0.0,
        'alpha': 4.73884,
        'beta': 0.15
    },
    'rotating_plasma_cylinder': {
        'k2': 1.0,
        'k3': 0.0,
        'cte_p0': 0.1,
        'p1': 8.0,
        'p2': 0.0,
        'p3': 0.0,
        'p4': 1.0,
        'p5': 0.0,
        'p6': 0.1
    },
    'rotating_theta_pinch': {
        'k2': 1.0,
        'k3': 0.0,
        'cte_rho0': 1.0,
        'alpha': 2.0,
        'p1': 0.1667,
        'p2': 1.0,
        'p3': 0.0
    },
    'suydam_cluster': {
        'k2': 1.0,
        'k3': -1.2,
        'alpha': 2.0,
        'cte_p0': 0.05,
        'cte_v03': 0.14,
        'p1': 0.1
    },
    'uniform_conduction': {
        'k2': 1.0,
        'k3': np.pi,
        'beta': 0.25
    }
}

PRECODED_MULTIRUNS = {
    'gravito_acoustic_NB10': {
        'equilibrium': 'gravito_acoustic',
        'gridpts': 31,
        'nb_runs': 40,
        'parameters': {
            'k2': np.linspace(0, np.sqrt(250), 40),
            'k3': np.linspace(0, np.sqrt(250), 40)
        }
    },
    'gravito_mhd_beta1': {
        'equilibrium': 'gravito_mhd',
        'gridpts': 31,
        'nb_runs': 40,
        'parameters': {
            'k2': np.linspace(0, np.sqrt(250), 40),
            'k3': np.linspace(0, np.sqrt(250), 40),
            'cte_p0': 0.5,  # so beta = 1
            'alpha': 20.0
        },
    },
    'gravito_mhd_beta50': {
        'equilibrium': 'gravito_mhd',
        'gridpts': 31,
        'nb_runs': 40,
        'parameters': {
            'k2': np.linspace(0, np.sqrt(250), 40),
            'k3': np.linspace(0, np.sqrt(250), 40),
            'cte_p0': 25.0, # so beta = 50
            'alpha': 20.0
        },
    },
    'interchange_modes_noshear': {
        'equilibrium': 'interchange_modes',
        'gridpts': 31,
        'nb_runs': 40,
        'parameters': {
            'k2': np.pi * np.sin(np.linspace(0, np.pi, 40)),
            'k3': np.pi * np.cos(np.linspace(0, np.pi, 40)),
            'lambda': 0.0,
            'cte_p0': 0.25, # so beta = 0.5
            'alpha': 20.0
        }
    },
    'interchange_modes_0.3shear': {
        'equilibrium': 'interchange_modes',
        'gridpts': 31,
        'nb_runs': 40,
        'parameters': {
            'k2': np.pi * np.sin(np.linspace(0, np.pi, 40)),
            'k3': np.pi * np.cos(np.linspace(0, np.pi, 40)),
            'lambda': 0.3,
            'cte_p0': 0.25, # so beta = 0.5
            'alpha': 20.0
        }
    },
    'constant_current_m-2': {
        'equilibrium': 'constant_current_tokamak',
        'gridpts': 51,
        'nb_runs': 39,
        'parameters': {
            'k2': -2.0,
            'k3': 0.2,
            'j0': (2.0 * 0.2) / np.linspace(1.9, 2.1, 39)
        }
    }
}
