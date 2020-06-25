import pytest
import pylbo
from pathlib import Path
from .suite_utils import get_filepaths, \
    get_answer_filepaths, \
    output, \
    FP_LIMIT

datfile, logfile = get_filepaths('resistive_homo')
answer_datfile, answer_logfile = get_answer_filepaths('resistive_homo')

config = {
    'geometry': 'Cartesian',
    'x_start': 0,
    'x_end': 1,
    'gridpoints': 51,
    'parameters': {
        'k2': 0.0,
        'k3': 1.0,
        'beta': 0.25
    },
    'equilibrium_type': 'resistive_homo',
    'resistivity': True,
    'use_fixed_resistivity': True,
    'fixed_eta_value': 0.001,
    'logging_level': 0,
    'show_results': False,
    'write_eigenfunctions': True,
    'write_matrices': False,
    'basename_datfile': datfile.stem,
    'basename_logfile': logfile.stem,
    'output_folder': str(output)
}

pylbo.set_loglevel('warning')

ev_guesses = [-0.4180167-0.0008883j,
              -0.4159505-0.0034422j,
              -0.4154890-0.0076967j,
              -0.4151935-0.0136528j,
              -0.4148061-0.0213107]


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


@pytest.fixture(scope='module', autouse=True)
def efs_test(ds_test):
    idxs, _ = ds_test.get_nearest_eigenvalues(ev_guesses)
    efh = pylbo.EigenfunctionHandler(ds_test)
    eigenfuncs = efh.get_eigenfunctions(idxs)
    return eigenfuncs


@pytest.fixture(scope='module', autouse=True)
def efs_answer(ds_answer):
    idxs, _ = ds_answer.get_nearest_eigenvalues(ev_guesses)
    efh = pylbo.EigenfunctionHandler(ds_answer)
    eigenfuncs = efh.get_eigenfunctions(idxs)
    return eigenfuncs


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


def test_rho_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        assert ef_test.get('rho') == pytest.approx(ef_answer.get('rho'), 1e-6)


def test_v1_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        assert ef_test.get('v1') == pytest.approx(ef_answer.get('v1'), 1e-6)
        # v1 nodes should be zero for wall conditions on edges
        for edge in (0, -1):
            assert ef_test.get('v1')[edge] == pytest.approx(0, FP_LIMIT)


def test_v2_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        assert ef_test.get('v2') == pytest.approx(ef_answer.get('v2'), 1e-6)


def test_v3_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        assert ef_test.get('v3') == pytest.approx(ef_answer.get('v3'), 1e-6)


def test_T_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        assert ef_test.get('T') == pytest.approx(ef_answer.get('T'), 1e-6)


def test_a1_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        assert ef_test.get('a1') == pytest.approx(ef_answer.get('a1'), 1e-6)


def test_a2_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        assert ef_test.get('a2') == pytest.approx(ef_answer.get('a2'), 1e-6)
        # v1 nodes should be zero for wall conditions on edges
        for edge in (0, -1):
            assert ef_test.get('a2')[edge] == pytest.approx(0, FP_LIMIT)


def test_a3_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        assert ef_test.get('a3') == pytest.approx(ef_answer.get('a3'), 1e-6)
        # v1 nodes should be zero for wall conditions on edges
        for edge in (0, -1):
            assert ef_test.get('a3')[edge] == pytest.approx(0, FP_LIMIT)
