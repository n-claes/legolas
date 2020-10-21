import pytest
import pylbo
import copy
from pathlib import Path
from .suite_utils import (
    get_filepaths,
    get_answer_filepaths,
    output,
    compare_eigenvalues,
    compare_eigenfunctions,
)

name = "resistive_homo"
datfile, logfile = get_filepaths(name)
answer_datfile, answer_logfile = get_answer_filepaths(name)

config = {
    "geometry": "Cartesian",
    "x_start": 0,
    "x_end": 1,
    "gridpoints": 51,
    "parameters": {
        "k2": 0.0,
        "k3": 1.0,
        "beta": 0.25,
        "cte_rho0": 1.0,
        "cte_B02": 0.0,
        "cte_B03": 1.0,
    },
    "equilibrium_type": "resistive_homo",
    "resistivity": True,
    "use_fixed_resistivity": True,
    "fixed_eta_value": 0.001,
    "logging_level": 0,
    "show_results": False,
    "write_eigenfunctions": True,
    "write_matrices": False,
    "basename_datfile": datfile.stem,
    "basename_logfile": logfile.stem,
    "output_folder": str(output),
}

pylbo.set_loglevel("warning")
ev_guesses = [
    -0.4180167 - 0.0008883j,
    -0.4159505 - 0.0034422j,
    -0.4154890 - 0.0076967j,
    -0.4151935 - 0.0136528j,
    -0.4148061 - 0.0213107,
]


@pytest.fixture(scope="module", autouse=True)
def ds_test():
    parfile = pylbo.generate_parfiles(
        parfile_dict=config, basename_parfile=datfile.stem, output_dir=output
    )
    pylbo.run_legolas(parfile, remove_parfiles=True)
    return pylbo.load(datfile)


@pytest.fixture(scope="module", autouse=True)
def ds_answer():
    return pylbo.load(answer_datfile)


@pytest.fixture(scope="module", autouse=True)
def log_test():
    return pylbo.read_log_file(logfile, sort=True)


@pytest.fixture(scope="module", autouse=True)
def log_answer():
    return pylbo.read_log_file(answer_logfile, sort=True)


@pytest.fixture(scope="module", autouse=True)
def efs_test(ds_test):
    idxs, _ = ds_test.get_nearest_eigenvalues(ev_guesses)
    efh = pylbo.EigenfunctionHandler(ds_test)
    eigenfuncs = efh.get_eigenfunctions(idxs)
    return eigenfuncs


@pytest.fixture(scope="module", autouse=True)
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


def test_params(ds_test, ds_answer):
    params = copy.deepcopy(ds_test.parameters)
    answ_params = copy.deepcopy(ds_answer.parameters)
    assert params.pop("k2") == pytest.approx(0) == answ_params.pop("k2")
    assert params.pop("k3") == pytest.approx(1) == answ_params.pop("k3")
    assert params.pop("beta") == pytest.approx(0.25) == answ_params.pop("beta")
    assert params.pop("cte_rho0") == pytest.approx(1) == answ_params.pop("cte_rho0")
    assert params.pop("cte_B02") == pytest.approx(0) == answ_params.pop("cte_B02")
    assert params.pop("cte_B03") == pytest.approx(1) == answ_params.pop("cte_B03")
    assert len(params) == 0
    assert len(answ_params) == 0


def test_eq_type(ds_test):
    assert ds_test.eq_type == "resistive_homo"


def test_eta_value(ds_test):
    for eta_val in ds_test.equilibria.get("eta"):
        assert eta_val == pytest.approx(1e-3)


def test_compare_evs(log_test, log_answer):
    compare_eigenvalues(log_test, log_answer, ds_name=name)


def test_rho_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        compare_eigenfunctions(ef_test.get("rho"), ef_answer.get("rho"), use_abs=True)


def test_v1_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        compare_eigenfunctions(ef_test.get("v1"), ef_answer.get("v1"), use_abs=True)


def test_v1_eigenfunction_ends(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        # v1 nodes should be zero for wall conditions on edges
        for edge in (0, -1):
            assert ef_test.get("v1").real[edge] == pytest.approx(0)
            assert ef_test.get("v1").imag[edge] == pytest.approx(0)


def test_v2_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        compare_eigenfunctions(ef_test.get("v2"), ef_answer.get("v2"), use_abs=True)


def test_v3_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        compare_eigenfunctions(ef_test.get("v3"), ef_answer.get("v3"), use_abs=True)


def test_T_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        compare_eigenfunctions(ef_test.get("T"), ef_answer.get("T"), use_abs=True)


def test_a1_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        compare_eigenfunctions(ef_test.get("a1"), ef_answer.get("a1"), use_abs=True)


def test_a2_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        compare_eigenfunctions(ef_test.get("a2"), ef_answer.get("a2"), use_abs=True)


def test_a2_eigenfunction_ends(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        # a2 nodes should be zero for wall conditions on edges
        for edge in (0, -1):
            assert ef_test.get("a2").real[edge] == pytest.approx(0)
            assert ef_test.get("a2").imag[edge] == pytest.approx(0)


def test_a3_eigenfunctions(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        compare_eigenfunctions(ef_test.get("a3"), ef_answer.get("a3"), use_abs=True)


def test_a3_eigenfunction_ends(efs_test, efs_answer):
    for ef_test, ef_answer in zip(efs_test, efs_answer):
        # a3 nodes should be zero for wall conditions on edges
        for edge in (0, -1):
            assert ef_test.get("a3").real[edge] == pytest.approx(0)
            assert ef_test.get("a3").imag[edge] == pytest.approx(0)
