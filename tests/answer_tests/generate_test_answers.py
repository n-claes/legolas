import pylbo
import copy
import os
import shutil
from pathlib import Path
from test_adiabatic_homo import config as adiabatic_config

output_dir = (pylbo.LEGOLAS_DIR / 'tests/answer_tests/answers').resolve()

tests = [adiabatic_config]

def main():
    for test_dict in tests:
        config = copy.deepcopy(test_dict)
        print('generating {}'.format(config['equilibrium_type']))
        new_basename_datfile = '{}_{}'.format('answer', config['basename_datfile'])
        new_basename_logfile = '{}_{}'.format('answer', config['basename_logfile'])
        config.update({'basename_datfile': new_basename_datfile,
                       'basename_logfile': new_basename_logfile,
                       'output_folder': str(output_dir)})
        parfile = pylbo.generate_parfiles(parfile_dict=config, basename_parfile=new_basename_datfile,
                                          output_dir=output_dir)
        pylbo.run_legolas(parfile, remove_parfiles=True)

if __name__ == '__main__':
    main()
