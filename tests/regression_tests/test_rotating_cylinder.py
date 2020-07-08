import pytest
import pylbo
import copy
from pathlib import Path
from .suite_utils import get_filepaths, \
    get_answer_filepaths, \
    output, \
    compare_eigenvalues

name = 'rotating_cylinder'
datfile, logfile = get_filepaths(name)
answer_datfile, answer_logfile = get_answer_filepaths(name)

config = {
    'geometry': 'cylindrical',
    'x_start': 0,
    'x_end': 1,
    'gridpoints': 51,
    'parameters': {
        'k2': 1.0,
        'k3': 0.0,
        'p1': 8.0,
        'p2': 0.0,
        'p3': 0.0,
        'p4': 1.0,
        'p5': 0.0,
        'p6': 0.0,
        'cte_p0': 0.1,
        'cte_rho0': 1
    },
    'equilibrium_type': 'rotating_plasma_cylinder',
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
    assert params.pop('k2') == pytest.approx(1)
    assert params.pop('k3') == pytest.approx(0)
    assert params.pop('cte_p0') == pytest.approx(0.1)
    assert params.pop('cte_rho0') == pytest.approx(1)
    assert params.pop('p1') == pytest.approx(8)
    assert params.pop('p2') == pytest.approx(0)
    assert params.pop('p3') == pytest.approx(0)
    assert params.pop('p4') == pytest.approx(1)
    assert params.pop('p5') == pytest.approx(0)
    assert params.pop('p6') == pytest.approx(0)
    assert len(params) == 0


def test_eq_type(ds_test):
    assert ds_test.eq_type == 'rotating_plasma_cylinder'


def test_compare_evs(log_test, log_answer):
    compare_eigenvalues(log_test, log_answer, ds_name=name)
