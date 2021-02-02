import pytest
import numpy as np
from pathlib import Path
import pylbo

from regression_tests.suite_utils import (
    output,
    get_filepaths,
    get_answer_filepaths,
    compare_eigenvalues,
    compare_eigenfunctions,
)
from regression_tests.test_adiabatic_homo import adiabatic_homo_setup
from regression_tests.test_discrete_alfven import discrete_alfven_setup
from regression_tests.test_fluxtube_coronal import fluxtube_coronal_setup
from regression_tests.test_fluxtube_photospheric import fluxtube_photospheric_config
from regression_tests.test_gold_hoyle import gold_hoyle_setup
from regression_tests.test_interchange_modes import interchange_modes_setup
from regression_tests.test_internal_kink import internal_kink_setup
from regression_tests.test_kh_cd import kh_cd_setup
from regression_tests.test_KHI import khi_setup
from regression_tests.test_magnetothermal import magnetothermal_setup
from regression_tests.test_MRI import mri_setup
from regression_tests.test_quasimodes import quasimodes_setup
from regression_tests.test_resistive_homo import resistive_homo_setup
from regression_tests.test_resistive_tearing import resistive_tearing_setup
from regression_tests.test_resistive_tearing_flow import resistive_tearing_flow_setup
from regression_tests.test_rotating_cylinder import rotating_cylinder_setup
from regression_tests.test_RTI import rti_setup
from regression_tests.test_RTI_theta_pinch_HD import rti_thetapinch_hd_setup
from regression_tests.test_RTI_theta_pinch_MHD import rti_thetapinch_mhd_setup
from regression_tests.test_suydam import suydam_setup
from regression_tests.test_tokamak import tokamak_setup


pylbo.set_loglevel("warning")

# retrieve test setups
tests_to_run = [
    adiabatic_homo_setup,
    discrete_alfven_setup,
    fluxtube_coronal_setup,
    fluxtube_photospheric_config,
    gold_hoyle_setup,
    interchange_modes_setup,
    internal_kink_setup,
    kh_cd_setup,
    khi_setup,
    magnetothermal_setup,
    mri_setup,
    quasimodes_setup,
    resistive_homo_setup,
    resistive_tearing_setup,
    resistive_tearing_flow_setup,
    rotating_cylinder_setup,
    rti_setup,
    rti_thetapinch_hd_setup,
    rti_thetapinch_mhd_setup,
    suydam_setup,
    tokamak_setup
]
# configure test setup
for _setup in tests_to_run:
    # get filepaths
    name = _setup["name"]
    datfile, logfile = get_filepaths(name)
    answer_datfile, answer_logfile = get_answer_filepaths(name)
    # set filepaths
    _setup["datfile"] = datfile
    _setup["logfile"] = logfile
    _setup["answer_datfile"] = answer_datfile
    _setup["answer_logfile"] = answer_logfile
    _setup["config"]["basename_datfile"] = datfile.stem
    _setup["config"]["basename_logfile"] = logfile.stem
    _setup["config"]["output_folder"] = str(output)
    # avoid running the same test multipe times
    _setup["test_needs_run"] = True
# set test IDs
ids = [_setup["name"] for _setup in tests_to_run]


# ===== GENERAL TESTS =====
@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_generation(ds_test, ds_answer, setup):
    assert ds_test.datfile == setup["datfile"]
    assert ds_answer.datfile == setup["answer_datfile"]


@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_existence(setup):
    assert Path.is_file(setup["datfile"])
    assert Path.is_file(setup["logfile"])
    assert Path.is_file(setup["answer_datfile"])
    assert Path.is_file(setup["answer_logfile"])


@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_eq_type(ds_test, setup):
    assert ds_test.eq_type == setup["config"]["equilibrium_type"]


@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_parameters(ds_test, ds_answer, setup):
    assert ds_test.parameters == ds_answer.parameters


@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_eigenvalues(log_test, log_answer, setup):
    compare_eigenvalues(log_test, log_answer, ds_name=setup["name"])


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("all_eigenvalues_real", False)],
    ids=[s["name"] for s in tests_to_run if s.get("all_eigenvalues_real", False)],
)
def test_if_eigenvalues_all_real(ds_test, setup):
    assert np.all(ds_test.eigenvalues.imag == pytest.approx(0))


# ===== TESTS FOR EIGENFUNCTIONS =====
@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_rho_eigenfunction(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        compare_eigenfunctions(ef_test.get("rho"), ef_answer.get("rho"), use_abs=True)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_v1_eigenfunction(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        compare_eigenfunctions(ef_test.get("v1"), ef_answer.get("v1"), use_abs=True)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_v1_eigenfunction_edges(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        # v1 must be zero on edges for wall boundary conditions
        for edge in (0, -1):
            assert ef_test.get("v1").real[edge] == pytest.approx(0)
            assert ef_test.get("v1").imag[edge] == pytest.approx(0)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_v2_eigenfunction(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        compare_eigenfunctions(ef_test.get("v2"), ef_answer.get("v2"), use_abs=True)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_v3_eigenfunction(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        compare_eigenfunctions(ef_test.get("v3"), ef_answer.get("v3"), use_abs=True)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_T_eigenfunction(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        compare_eigenfunctions(ef_test.get("T"), ef_answer.get("T"), use_abs=True)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_a1_eigenfunction(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        compare_eigenfunctions(ef_test.get("a1"), ef_answer.get("a1"), use_abs=True)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_a2_eigenfunction(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        compare_eigenfunctions(ef_test.get("a2"), ef_answer.get("a2"), use_abs=True)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_a2_eigenfunction_edges(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        # a2 must be zero on edges for wall boundary conditions
        for edge in (0, -1):
            assert ef_test.get("a2").real[edge] == pytest.approx(0)
            assert ef_test.get("a2").imag[edge] == pytest.approx(0)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_a3_eigenfunction(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        compare_eigenfunctions(ef_test.get("a3"), ef_answer.get("a3"), use_abs=True)


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("ev_guesses") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("ev_guesses") is not None],
)
def test_a3_eigenfunction_edges(eigfuncs_test, eigfuncs_answer, setup):
    for ef_test, ef_answer in zip(eigfuncs_test, eigfuncs_answer):
        # a2 must be zero on edges for wall boundary conditions
        for edge in (0, -1):
            assert ef_test.get("a3").real[edge] == pytest.approx(0)
            assert ef_test.get("a3").imag[edge] == pytest.approx(0)
