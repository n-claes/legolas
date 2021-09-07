import pytest
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images
from pathlib import Path
import pylbo

from regression_tests.suite_utils import (
    output,
    get_filepaths,
    get_answer_filepaths,
    get_image_filename,
    ABS_TOL,
    SAVEFIG_KWARGS,
    RMS_TOLERANCE,
    EF_NAMES,
)
from regression_tests.test_adiabatic_homo import adiabatic_homo_setup
from regression_tests.test_adiabatic_homo_ef_subset import adiabatic_ef_subset_setup
from regression_tests.test_couette_flow import couette_flow_setup
from regression_tests.test_discrete_alfven import discrete_alfven_setup
from regression_tests.test_fluxtube_coronal import fluxtube_coronal_setup
from regression_tests.test_fluxtube_photospheric import fluxtube_photospheric_setup
from regression_tests.test_gold_hoyle import gold_hoyle_setup
from regression_tests.test_interchange_modes import interchange_modes_setup
from regression_tests.test_internal_kink import internal_kink_setup
from regression_tests.test_kh_cd import kh_cd_setup
from regression_tests.test_KHI import khi_setup
from regression_tests.test_magnetothermal import magnetothermal_setup
from regression_tests.test_magnetothermal_arnoldi_si import magneto_arnoldi_si_setup
from regression_tests.test_MRI import mri_setup
from regression_tests.test_quasimodes import quasimodes_setup
from regression_tests.test_resistive_homo import resistive_homo_setup
from regression_tests.test_resistive_homo_arnoldi_si import (
    resistive_homo_arnoldi_si_setup,
)
from regression_tests.test_resistive_homo_ef_subset import (
    resistive_homo_ef_subset_setup,
)
from regression_tests.test_resistive_tearing import resistive_tearing_setup
from regression_tests.test_resistive_tearing_flow import resistive_tearing_flow_setup
from regression_tests.test_rotating_cylinder import rotating_cylinder_setup
from regression_tests.test_RTI import rti_setup
from regression_tests.test_RTI_KHI import rti_khi_setup
from regression_tests.test_RTI_theta_pinch_HD import rti_thetapinch_hd_setup
from regression_tests.test_RTI_theta_pinch_MHD import rti_thetapinch_mhd_setup
from regression_tests.test_suydam import suydam_setup
from regression_tests.test_taylor_couette import taylor_couette_setup
from regression_tests.test_tokamak import tokamak_setup

from regression_tests.test_multi_constant_current import constant_current_setup
from regression_tests.test_multi_coronal_fluxtube import coronal_tube_setup
from regression_tests.test_multi_gravito_HD import gravito_hd_setup
from regression_tests.test_multi_gravito_MHD import gravito_mhd_setup
from regression_tests.test_multi_interchange import interchange_setup
from regression_tests.test_multi_iso_atmo import iso_atmo_beta_half_setup
from regression_tests.test_multi_photospheric_fluxtube import photospheric_tube_setup


# retrieve test setups
tests_to_run = [
    adiabatic_homo_setup,
    adiabatic_ef_subset_setup,
    couette_flow_setup,
    discrete_alfven_setup,
    fluxtube_coronal_setup,
    fluxtube_photospheric_setup,
    gold_hoyle_setup,
    interchange_modes_setup,
    internal_kink_setup,
    kh_cd_setup,
    khi_setup,
    magnetothermal_setup,
    magneto_arnoldi_si_setup,
    mri_setup,
    quasimodes_setup,
    resistive_homo_setup,
    resistive_homo_arnoldi_si_setup,
    resistive_homo_ef_subset_setup,
    resistive_tearing_setup,
    resistive_tearing_flow_setup,
    rotating_cylinder_setup,
    rti_setup,
    rti_khi_setup,
    rti_thetapinch_hd_setup,
    rti_thetapinch_mhd_setup,
    suydam_setup,
    taylor_couette_setup,
    tokamak_setup,
]
multirun_tests = [
    constant_current_setup,
    coronal_tube_setup,
    gravito_hd_setup,
    gravito_mhd_setup,
    interchange_setup,
    iso_atmo_beta_half_setup,
    photospheric_tube_setup,
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
# configure multirun test setup
for _setup in multirun_tests:
    name = _setup["name"]
    answer_series, _ = get_answer_filepaths(name)
    _setup["answer_series"] = answer_series.with_suffix(".pickle")
    _setup["test_needs_run"] = True
    _setup["config"]["basename_datfile"] = f"{name}_series"
    _setup["config"]["output_folder"] = str(output)

# set test IDs
ids = [_setup["name"] for _setup in tests_to_run]
multirun_ids = [_setup["name"] for _setup in multirun_tests]


# ===== EXISTENCE OF ANSWER FILES =====
# single spectra
@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_answer_datfile_exists(setup):
    if not setup["answer_datfile"].is_file():
        raise FileNotFoundError(f"{setup['answer_datfile']} does not exist!")


# multirun spectra
@pytest.mark.parametrize("setup", multirun_tests, ids=multirun_ids)
def test_answer_multirun_file_exists(setup):
    if not setup["answer_series"].is_file():
        raise FileNotFoundError(f"{setup['answer_series']} does not exist!")


# ===== GENERATION OF SINGLE FILES =====
@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_generate_datfile(ds_test, ds_answer, setup):
    assert ds_test.datfile == setup["datfile"]
    assert ds_answer.datfile == setup["answer_datfile"]


# ===== EXISTENCE OF SINGLE FILES =====
@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_datfile_exists(setup):
    if not setup["datfile"].is_file():
        raise FileNotFoundError(f"{setup['datfile']} does not exist!")


# ===== GENERATION OF MULTIRUN FILES =====
@pytest.mark.parametrize("setup", multirun_tests, ids=multirun_ids)
def test_generate_multirun_file(series_test, series_answer, setup):
    assert len(series_test) == setup["config"]["number_of_runs"]
    assert len(series_test) == len(series_answer)


# ===== EXISTENCE OF MULTIRUN FILES =====
@pytest.mark.parametrize("setup", multirun_tests, ids=multirun_ids)
def test_multirun_file_exists(setup):
    files = sorted(
        Path(setup["config"]["output_folder"]).glob(
            f"*{setup['config']['basename_datfile']}.dat"
        )
    )
    if not files:
        raise FileNotFoundError(f"No multirun datfiles found for {setup['name']}!")


# ===== GENERATION OF SPECTRUM IMAGES =====
@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("image_limits") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("image_limits") is not None],
)
def test_generate_spectrum_images(ds_test, ds_answer, setup, imagedir):
    p_test = pylbo.plot_spectrum(ds_test)
    p_answer = pylbo.plot_spectrum(ds_answer)
    setup["spectrum_images"] = []
    for image_lims in setup.get("image_limits"):
        figname = f"{get_image_filename(setup['name'], limits_dict=image_lims)}.png"
        figname_answer = f"{figname[:-4]}-baseline.png"
        xlims = image_lims.get("xlims")
        ylims = image_lims.get("ylims")
        # save test image
        p_test.ax.set_xlim(xlims)
        p_test.ax.set_ylim(ylims)
        p_test.fig.savefig(imagedir / figname, **SAVEFIG_KWARGS)
        # save baseline image
        p_answer.ax.set_xlim(xlims)
        p_answer.ax.set_ylim(ylims)
        p_answer.fig.savefig(imagedir / figname_answer, **SAVEFIG_KWARGS)
        setup["spectrum_images"].append((figname, figname_answer))
    plt.close(p_test.fig)
    plt.close(p_answer.fig)


# ===== GENERATION OF EIGENFUNCTION IMAGES =====
@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("eigenfunctions") is not None],
    ids=[s["name"] for s in tests_to_run if s.get("eigenfunctions") is not None],
)
def test_generate_eigenfunction_images(ds_test, ds_answer, setup, imagedir):
    setup["eigenfunction_images"] = []
    fig_test, ax_test = plt.subplots(3, 3, figsize=(10, 10), sharex="all")
    axes_test = ax_test.flatten()
    fig_test.delaxes(axes_test[-1])
    fig_answer, ax_answer = plt.subplots(3, 3, figsize=(10, 10), sharex="all")
    axes_answer = ax_answer.flatten()
    fig_answer.delaxes(axes_answer[-1])
    for i, efs in enumerate(setup.get("eigenfunctions")):
        eigenval = efs.get("eigenvalue")
        figname = f"{setup['name']}_eigenfunctions_{i}.png"
        figtitle = f"{setup['name']} -- eigenvalue={eigenval:.6f}"
        # generate baseline image
        (efs_answer,) = ds_answer.get_eigenfunctions(ev_guesses=eigenval)
        figname_answer = f"{figname[:-4]}-baseline.png"
        for ax, ef_name in zip(axes_answer, EF_NAMES):
            result = efs_answer[ef_name].real + efs_answer[ef_name].imag
            # small values to zero
            result[np.where(abs(result) < ABS_TOL)] = 0
            ax.plot(ds_answer.ef_grid, abs(result), lw=3)
            ax.set_yticks([])
            ax.set_title(ef_name)
        fig_answer.suptitle(figtitle)
        fig_answer.tight_layout()
        fig_answer.savefig(imagedir / figname_answer, **SAVEFIG_KWARGS)
        [ax.clear() for ax in axes_answer]

        # generate test image
        (efs_test,) = ds_test.get_eigenfunctions(ev_guesses=eigenval)
        for ax, ef_name in zip(axes_test, EF_NAMES):
            result = efs_test[ef_name].real + efs_test[ef_name].imag
            # small values to zero
            result[np.where(abs(result) < ABS_TOL)] = 0
            ax.plot(ds_test.ef_grid, abs(result), lw=3)
            ax.set_yticks([])
            ax.set_title(ef_name)
        fig_test.suptitle(figtitle)
        fig_test.tight_layout()
        fig_test.savefig(imagedir / figname, **SAVEFIG_KWARGS)
        [ax.clear() for ax in axes_test]

        setup["eigenfunction_images"].append((figname, figname_answer))
    plt.close(fig_answer)
    plt.close(fig_test)


# ===== GENERATION OF DERIVED EIGENFUNCTION IMAGES =====
@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("derived_eigenfunctions") is not None],
    ids=[
        s["name"] for s in tests_to_run if s.get("derived_eigenfunctions") is not None
    ],
)
def test_generate_derived_eigenfunction_images(ds_test, ds_answer, setup, imagedir):
    setup["derived_eigenfunction_images"] = []
    fig_test, axes_test = plt.subplots(5, 4, figsize=(10, 8), sharex="all")
    axes_test = axes_test.flatten()
    fig_answer, axes_answer = plt.subplots(5, 4, figsize=(10, 8), sharex="all")
    axes_answer = axes_answer.flatten()
    for i, efs in enumerate(setup.get("derived_eigenfunctions")):
        eigenval = efs.get("eigenvalue")
        figname = f"{setup['name']}_derived_eigenfunctions_{i}.png"
        figtitle = f"{setup['name']} -- eigenvalue={eigenval:.6f}"
        # generate baseline image
        (efs_answer,) = ds_answer.get_derived_eigenfunctions(ev_guesses=eigenval)
        figname_answer = f"{figname[:-4]}-baseline.png"
        for ax, ef_name in zip(axes_answer, ds_answer.derived_ef_names):
            result = efs_answer[ef_name].real + efs_answer[ef_name].imag
            # small values to zero
            result[np.where(abs(result) < ABS_TOL)] = 0
            ax.plot(ds_answer.ef_grid, abs(result), lw=3)
            ax.set_yticks([])
            ax.set_title(ef_name)
        fig_answer.suptitle(figtitle)
        fig_answer.tight_layout()
        fig_answer.savefig(imagedir / figname_answer, **SAVEFIG_KWARGS)
        [ax.clear() for ax in axes_answer]

        # generate test image
        (efs_test,) = ds_test.get_derived_eigenfunctions(ev_guesses=eigenval)
        for ax, ef_name in zip(axes_test, ds_answer.derived_ef_names):
            result = efs_test[ef_name].real + efs_test[ef_name].imag
            # small values to zero
            result[np.where(abs(result) < ABS_TOL)] = 0
            ax.plot(ds_test.ef_grid, abs(result), lw=3)
            ax.set_yticks([])
            ax.set_title(ef_name)
        fig_test.suptitle(figtitle)
        fig_test.tight_layout()
        fig_test.savefig(imagedir / figname, **SAVEFIG_KWARGS)
        [ax.clear() for ax in axes_test]

        setup["derived_eigenfunction_images"].append((figname, figname_answer))
    plt.close(fig_answer)
    plt.close(fig_test)


# ===== GENERATION OF MULTIRUN IMAGES =====
@pytest.mark.parametrize("setup", multirun_tests, ids=multirun_ids)
def test_generate_multi_spectra_images(series_test, series_answer, setup, imagedir):
    setup["spectrum_images"] = []
    p_test = pylbo.plot_spectrum_multi(
        series_test, xdata=setup["xdata"], use_squared_omega=setup["use_squared_omega"]
    )
    p_answer = pylbo.plot_spectrum_multi(
        series_answer,
        xdata=setup["xdata"],
        use_squared_omega=setup["use_squared_omega"],
    )
    xlims = setup["limits"]["xlims"]
    ylims = setup["limits"]["ylims"]
    for p in (p_test, p_answer):
        p.set_x_scaling(setup.get("x_scaling", 1))
        p.set_y_scaling(setup.get("y_scaling", 1))
        if setup.get("symlog", None) is not None:
            p.ax.set_yscale("symlog", linthresh=setup["symlog"])
        p.ax.set_xlim(xlims)
        p.ax.set_ylim(ylims)
    figname = f"{get_image_filename(setup['name'], limits_dict=setup['limits'])}.png"
    figname_answer = f"{figname[:-4]}-baseline.png"
    # save test image
    p_test.fig.savefig(imagedir / figname, **SAVEFIG_KWARGS)
    # save baseline image
    p_answer.fig.savefig(imagedir / figname_answer, **SAVEFIG_KWARGS)
    setup["spectrum_images"].append((figname, figname_answer))
    plt.close(p_test.fig)
    plt.close(p_answer.fig)


@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_eq_type(ds_test, setup):
    assert ds_test.eq_type == setup["config"]["equilibrium_type"]


@pytest.mark.parametrize("setup", tests_to_run, ids=ids)
def test_parameters(ds_test, ds_answer, setup):
    assert ds_test.parameters == ds_answer.parameters


# ===== TESTS FOR EIGENVALUES AND SPECTRUM =====
@pytest.mark.parametrize(
    argnames="setup,idx",
    argvalues=[
        (_s, _idx)
        for _s in tests_to_run
        if _s.get("image_limits") is not None
        for _idx, _ in enumerate(_s["image_limits"])
    ],
    ids=[
        f"{_s['name']}-[x={lims['xlims']}, y={lims['ylims']}]"
        for _s in tests_to_run
        if _s.get("image_limits") is not None
        for _idx, lims in enumerate(_s["image_limits"])
    ],
)
def test_eigenvalue_spectrum(imagedir, setup, idx, keep_files):
    test_image, baseline_image = setup["spectrum_images"][idx]
    result = compare_images(
        str(imagedir / baseline_image),
        str(imagedir / test_image),
        tol=setup["image_limits"][idx].get("RMS_TOLERANCE", RMS_TOLERANCE),
    )
    # result will be None if test succeeds, if pass we remove the images
    if result is not None:
        pytest.fail(result, pytrace=False)
    else:
        if not keep_files:
            Path(imagedir / baseline_image).unlink()
            Path(imagedir / test_image).unlink()


@pytest.mark.parametrize("setup", multirun_tests, ids=multirun_ids)
def test_multirun_spectrum(imagedir, setup, keep_files):
    test_image, baseline_image = setup["spectrum_images"][0]
    result = compare_images(
        str(imagedir / baseline_image),
        str(imagedir / test_image),
        tol=setup.get("RMS_TOLERANCE", RMS_TOLERANCE),
    )
    if result is not None:
        pytest.fail(result, pytrace=False)
    else:
        if not keep_files:
            Path(imagedir / baseline_image).unlink()
            Path(imagedir / test_image).unlink()


@pytest.mark.parametrize(
    "setup",
    [s for s in tests_to_run if s.get("all_eigenvalues_real", False)],
    ids=[s["name"] for s in tests_to_run if s.get("all_eigenvalues_real", False)],
)
def test_if_eigenvalues_all_real(ds_test, setup):
    assert np.all(ds_test.eigenvalues.imag == pytest.approx(0))


# ===== TESTS FOR EIGENFUNCTIONS =====
@pytest.mark.parametrize(
    argnames="setup,idx",
    argvalues=[
        (_s, _idx)
        for _s in tests_to_run
        if _s.get("eigenfunctions") is not None
        for _idx, _ in enumerate(_s["eigenfunctions"])
    ],
    ids=[
        f"{_s['name']}-eigenvalue={efs['eigenvalue']}"
        for _s in tests_to_run
        if _s.get("eigenfunctions") is not None
        for _idx, efs in enumerate(_s["eigenfunctions"])
    ],
)
def test_eigenfunction(imagedir, setup, idx, keep_files):
    test_image, baseline_image = setup["eigenfunction_images"][idx]
    result = compare_images(
        str(imagedir / baseline_image),
        str(imagedir / test_image),
        tol=setup["eigenfunctions"][idx].get("RMS_TOLERANCE", RMS_TOLERANCE),
    )
    # result will be None if test succeeds, if pass we remove images
    if result is not None:
        pytest.fail(result, pytrace=False)
    else:
        if not keep_files:
            Path(imagedir / baseline_image).unlink()
            Path(imagedir / test_image).unlink()


@pytest.mark.parametrize(
    argnames="setup,idx",
    argvalues=[
        (_s, _idx)
        for _s in tests_to_run
        if _s.get("eigenfunctions") is not None
        for _idx, _ in enumerate(_s["eigenfunctions"])
    ],
    ids=[
        f"{_s['name']}-eigenvalue={efs['eigenvalue']}"
        for _s in tests_to_run
        if _s.get("eigenfunctions") is not None
        for _idx, efs in enumerate(_s["eigenfunctions"])
    ],
)
def test_eigenfunction_edges(ds_test, ds_answer, setup, idx):
    eigenvalue = setup.get("eigenfunctions")[idx].get("eigenvalue")
    (efs_test,) = ds_test.get_eigenfunctions(ev_guesses=eigenvalue)
    (efs_answer,) = ds_answer.get_eigenfunctions(ev_guesses=eigenvalue)
    # v1, a2 and a3 must be zero on edges for wall boundary conditions
    for ef_name in ("v1", "a2", "a3"):
        for edge in (0, -1):
            assert efs_test.get(ef_name).real[edge] == pytest.approx(0)
            assert efs_test.get(ef_name).imag[edge] == pytest.approx(0)
            assert efs_answer.get(ef_name).real[edge] == pytest.approx(0)
            assert efs_answer.get(ef_name).imag[edge] == pytest.approx(0)


# ===== TESTS FOR DERIVED EIGENFUNCTION QUANTITIES =====
@pytest.mark.parametrize(
    argnames="setup,idx",
    argvalues=[
        (_s, _idx)
        for _s in tests_to_run
        if _s.get("derived_eigenfunctions") is not None
        for _idx, _ in enumerate(_s["derived_eigenfunctions"])
    ],
    ids=[
        f"{_s['name']}-eigenvalue={efs['eigenvalue']}"
        for _s in tests_to_run
        if _s.get("derived_eigenfunctions") is not None
        for _idx, efs in enumerate(_s["derived_eigenfunctions"])
    ],
)
def test_derived_eigenfunction(imagedir, setup, idx, keep_files):
    test_image, baseline_image = setup["derived_eigenfunction_images"][idx]
    result = compare_images(
        str(imagedir / baseline_image),
        str(imagedir / test_image),
        tol=setup["derived_eigenfunctions"][idx].get("RMS_TOLERANCE", RMS_TOLERANCE),
    )
    # result will be None if test succeeds, if pass we remove images
    if result is not None:
        pytest.fail(result, pytrace=False)
    else:
        if not keep_files:
            Path(imagedir / baseline_image).unlink()
            Path(imagedir / test_image).unlink()
