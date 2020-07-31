import pytest
import pylbo
import copy
from pathlib import Path
from .suite_utils import get_filepaths, \
    get_answer_filepaths, \
    output, \
    compare_eigenvalues

name = 'RTI_theta_pinch_MHD'
datfile, logfile = get_filepaths(name)
answer_datfile, answer_logfile = get_answer_filepaths(name)

config = {
    'gridpoints': 51,
    'parameters': {
        'k2': 1.0,
        'k3': 0.1,
        'cte_rho0': 1.0,
        'alpha': 2.0,
        'delta': 1/6,
        'r0': 0.0
    },
    'flow': True,
    'equilibrium_type': 'RTI_theta_pinch',
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
    assert params.pop('k2') == pytest.approx(1.0) == answ_params.pop('k2')
    assert params.pop('k3') == pytest.approx(0.1) == answ_params.pop('k3')
    assert params.pop('cte_rho0') == pytest.approx(1.0) == answ_params.pop('cte_rho0')
    assert params.pop('cte_p0') == pytest.approx(25/72) == answ_params.pop('cte_p0')
    assert params.pop('alpha') == pytest.approx(2.0) == answ_params.pop('alpha')
    assert params.pop('delta') == pytest.approx(1/6) == answ_params.pop('delta')
    assert params.pop('r0') == pytest.approx(0.0) == answ_params.pop('r0')
    assert len(params) == 0
    assert len(answ_params) == 0


def test_eq_type(ds_test):
    assert ds_test.eq_type == 'RTI_theta_pinch'


def test_compare_evs(log_test, log_answer):
    compare_eigenvalues(log_test, log_answer, ds_name=name)
