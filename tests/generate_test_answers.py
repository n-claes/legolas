import pylbo
import copy
from pathlib import Path
from regression_tests.test_adiabatic_homo import config as adiabatic_config
from regression_tests.test_discrete_alfven import config as discrete_alfven_config
from regression_tests.test_gold_hoyle import config as gold_hoyle_config
from regression_tests.test_interchange_modes import config as interchange_config
from regression_tests.test_internal_kink import config as internal_kink_config
from regression_tests.test_kh_cd import config as kh_cd_config
from regression_tests.test_KHI import config as KHI_config
from regression_tests.test_magnetothermal import config as magnetothermal_config
from regression_tests.test_quasimodes import config as quasi_config
from regression_tests.test_resistive_homo import config as resistive_homo_config
from regression_tests.test_resistive_tearing import config as resistive_tearing_config
from regression_tests.test_resistive_tearing_flow import config as resistive_tearing_flow_config
from regression_tests.test_rotating_cylinder import config as plasma_cylinder_config
from regression_tests.test_RTI import config as RTI_config
from regression_tests.test_RTI_KHI import config as RTI_KHI_config
from regression_tests.test_suydam import config as suydam_config
from regression_tests.test_tokamak import config as tokamak_config

pylbo.set_loglevel('info')
output_dir = (pylbo.LEGOLAS_DIR / 'tests/regression_tests/answers').resolve()

tests = [
    adiabatic_config,
    discrete_alfven_config,
    gold_hoyle_config,
    interchange_config,
    internal_kink_config,
    kh_cd_config,
    KHI_config,
    magnetothermal_config,
    quasi_config,
    resistive_homo_config,
    resistive_tearing_config,
    resistive_tearing_flow_config,
    plasma_cylinder_config,
    RTI_config,
    RTI_KHI_config,
    suydam_config,
    tokamak_config
]


def overwrite_files(base_filename):
    datfile = (output_dir / (base_filename + '.dat')).resolve()
    if Path.is_file(datfile):
        print('{} already exists!'.format(datfile.name))
        force = input('overwrite? ')
        if force.lower() in ('yes', 'y'):
            return True
        else:
            return False
    return True


def main():
    parfiles = []
    for test_dict in tests:
        config = copy.deepcopy(test_dict)
        name = config['equilibrium_type']
        print('='*50)
        print('>> generating {}'.format(name))
        new_basename_datfile = '{}_{}'.format('answer', config['basename_datfile'])
        new_basename_logfile = '{}_{}'.format('answer', config['basename_logfile'])
        config.update({'basename_datfile': new_basename_datfile,
                       'basename_logfile': new_basename_logfile,
                       'output_folder': str(output_dir)})
        if overwrite_files(new_basename_datfile):
            parfile = pylbo.generate_parfiles(parfile_dict=config,
                                              basename_parfile=new_basename_datfile,
                                              output_dir=output_dir)
            parfiles.append(*parfile)
        else:
            print('Skipping {}'.format(name))
    if parfiles:
        pylbo.run_legolas(parfiles, remove_parfiles=True, nb_cpus=4)


if __name__ == '__main__':
    main()
