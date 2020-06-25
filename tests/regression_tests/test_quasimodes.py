import pytest
import pylbo
from pathlib import Path
from .suite_utils import get_filepaths, \
    get_answer_filepaths, \
    output

datfile, logfile = get_filepaths('quasimodes')
answer_datfile, answer_logfile = get_answer_filepaths('quasimodes')

config = {
    'geometry': 'Cartesian',
    'x_start': 0.0,
    'x_end': 1.0,
    'gridpoints': 51,
    'parameters': {
        'k2': 1.0,
        'k3': 0.05,
        'p1': 0.9,
        'p2': 0.1,
        'r0': 0.1
    },
    'equilibrium_type': 'resonant_absorption',
    'resistivity': True,
    'use_fixed_resistivity': True,
    'fixed_eta_value': 10**(-4.0),
    'logging_level': 0,
    'show_results': False,
    'write_eigenfunctions': False,
    'write_matrices': False,
    'basename_datfile': datfile.stem,
    'basename_logfile': logfile.stem,
    'output_folder': str(output)
}

pylbo.set_loglevel('warning')


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
