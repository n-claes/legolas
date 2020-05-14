import pytest
import pylbo
from pathlib import Path

output = (pylbo.LEGOLAS_DIR / 'tests/answer_tests').resolve()
datfile = (output / 'adiabatic_homo.dat').resolve()
logfile = (output / 'adiabatic_homo.log').resolve()

answer_datfile = (output / 'answers/answer_{}'.format(datfile.name)).resolve()
answer_logfile = (output / 'answers/answer_{}'.format(logfile.name)).resolve()

config = {
    'gridpoints': 51,
    'geometry': 'Cartesian',
    'x_start': 0,
    'x_end': 1,
    'equilibrium_type': 'adiabatic_homo',
    'logging_level': 0,
    'show_results': False,
    'write_eigenfunctions': False,
    'write_matrices': False,
    'basename_datfile': datfile.stem,
    'basename_logfile': logfile.stem,
    'output_folder': str(output)
}

@pytest.fixture(scope='module', autouse=True)
def ds_test():
    parfile = pylbo.generate_parfiles(parfile_dict=config, basename_parfile=datfile.stem, output_dir=output)
    pylbo.run_legolas(parfile, remove_parfiles=True)
    return pylbo.load(datfile)

@pytest.fixture(scope='module', autouse=True)
def ds_answer():
    return pylbo.load(answer_datfile)

@pytest.fixture(scope='module', autouse=True)
def log_test():
    return pylbo.read_log_file(logfile, sort=True)

@pytest.fixture(scope='module', autouse=True)
def log_answer():
    return pylbo.read_log_file(answer_logfile, sort=True)

def test_existence():
    assert Path.is_file(datfile)
    assert Path.is_file(logfile)
    assert Path.is_file(answer_datfile)
    assert Path.is_file(answer_logfile)
    assert Path.is_dir(output)

def test_filenames(ds_test, ds_answer):
    assert ds_test.datfile == str(datfile)
    assert ds_answer.datfile == str(answer_datfile)

def test_compare_evs(log_test, log_answer):
    for i in range(len(log_test)):
        assert log_test[i] == log_answer[i]
