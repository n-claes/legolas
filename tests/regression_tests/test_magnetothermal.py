import pytest
import pylbo
import copy
from pathlib import Path
from .suite_utils import (
    get_filepaths,
    get_answer_filepaths,
    output,
    compare_eigenvalues,
)

name = "magnetothermal"
datfile, logfile = get_filepaths(name)
answer_datfile, answer_logfile = get_answer_filepaths(name)

config = {
    "geometry": "cylindrical",
    "x_start": 0,
    "x_end": 1,
    "gridpoints": 51,
    "parameters": {
        "k2": 0.0,
        "k3": 1.0,
        "cte_T0": 1.0,
    },
    "equilibrium_type": "magnetothermal_instabilities",
    "radiative_cooling": True,
    "cooling_curve": "rosner",
    "thermal_conduction": True,
    "use_fixed_tc_perp": True,
    "fixed_tc_perp_value": 0,
    "cgs_units": True,
    "unit_temperature": 2.6e6,
    "unit_magneticfield": 10.0,
    "unit_length": 1.0e8,
    "logging_level": 0,
    "show_results": False,
    "write_eigenfunctions": False,
    "write_matrices": False,
    "basename_datfile": datfile.stem,
    "basename_logfile": logfile.stem,
    "output_folder": str(output),
}

pylbo.set_loglevel("warning")


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
    assert params.pop("cte_T0") == pytest.approx(1) == answ_params.pop("cte_T0")
    assert len(params) == 0
    assert len(answ_params) == 0


def test_eq_type(ds_test):
    assert ds_test.eq_type == "magnetothermal_instabilities"


def test_conduction(ds_test):
    for kappa_perp_val in ds_test.equilibria.get("kappa_perp"):
        assert kappa_perp_val == pytest.approx(0)


def test_compare_evs(log_test, log_answer):
    compare_eigenvalues(log_test, log_answer, ds_name=name)


def test_units(ds_test):
    units = ds_test.units
    assert ds_test.cgs
    assert units.get("unit_temperature") == pytest.approx(2.6e6)
    assert units.get("unit_magneticfield") == pytest.approx(10)
    assert units.get("unit_length") == pytest.approx(1e8)
