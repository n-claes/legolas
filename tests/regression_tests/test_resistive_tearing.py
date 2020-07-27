import pytest
import pylbo
import copy
from pathlib import Path
from .suite_utils import get_filepaths, \
    get_answer_filepaths, \
    output, \
    compare_eigenvalues

name = 'resistive_tearing'
datfile, logfile = get_filepaths(name)
answer_datfile, answer_logfile = get_answer_filepaths(name)

config = {
    'geometry': 'Cartesian',
    'x_start': -0.5,
    'x_end': 0.5,
    'gridpoints': 51,
    'parameters': {
        'k2': 0.49,
        'k3': 0.0,
        'alpha': 4.73884,
        'beta': 0.15,
        'cte_rho0': 1.0
    },
    'equilibrium_type': 'resistive_tearing',
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


def test_params(ds_test, ds_answer):
    params = copy.deepcopy(ds_test.parameters)
    answ_params = copy.deepcopy(ds_answer.parameters)
    assert params.pop('k2') == pytest.approx(0.49) == answ_params.pop('k2')
    assert params.pop('k3') == pytest.approx(0) == answ_params.pop('k3')
    assert params.pop('alpha') == pytest.approx(4.73884) == answ_params.pop('alpha')
    assert params.pop('beta') == pytest.approx(0.15) == answ_params.pop('beta')
    assert params.pop('cte_rho0') == pytest.approx(1) == answ_params.pop('cte_rho0')
    assert len(params) == 0
    assert len(answ_params) == 0


def test_eq_type(ds_test):
    assert ds_test.eq_type == 'resistive_tearing'


def test_eta_value(ds_test):
    for eta_val in ds_test.equilibria.get('eta'):
        assert eta_val == pytest.approx(1e-4)


def test_compare_evs(log_test, log_answer):
    compare_eigenvalues(log_test, log_answer, ds_name=name)
