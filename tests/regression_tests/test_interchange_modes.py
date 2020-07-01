import pytest
import numpy as np
import pylbo
import copy
from pathlib import Path
from .suite_utils import get_filepaths, \
    get_answer_filepaths, \
    output, \
    compare_eigenvalues

name = 'interchange_modes'
datfile, logfile = get_filepaths(name)
answer_datfile, answer_logfile = get_answer_filepaths(name)

config = {
    'geometry': 'Cartesian',
    'x_start': 0,
    'x_end': 1,
    'gridpoints': 51,
    'parameters': {
        'k2': np.pi,
        'k3': np.pi,
        'g': 0.5,
        'cte_p0': 0.25,
        'lambda': 0,
        'alpha': 20.0
    },
    'equilibrium_type': 'interchange_modes',
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


def test_params(ds_test):
    params = copy.deepcopy(ds_test.parameters)
    assert params.pop('k2') == pytest.approx(np.pi)
    assert params.pop('k3') == pytest.approx(np.pi)
    assert params.pop('g') == pytest.approx(0.5)
    assert params.pop('cte_p0') == pytest.approx(0.25)
    assert params.pop('cte_rho0') == pytest.approx(30)
    assert params.pop('lambda') == pytest.approx(0)
    assert params.pop('alpha') == pytest.approx(20)
    assert params.pop('beta') == 0.5
    assert len(params) == 0


def test_eq_type(ds_test):
    assert ds_test.eq_type == 'interchange_modes'


def test_gravity_value(ds_test):
    for g_val in ds_test.equilibria.get('grav'):
        assert g_val == pytest.approx(0.5)


def test_compare_evs(log_test, log_answer):
    compare_eigenvalues(log_test, log_answer, ds_name=name)
