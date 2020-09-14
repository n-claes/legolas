import shutil
import numpy as np
import pytest
from pathlib import Path
import pylbo

pylbo.set_loglevel("warning")
FIG_DPI = 200
BASELINE_DIR = str(Path("utility_files/baseline").resolve())
RESULTS_DIR = str(Path("results").resolve())
NB_RUNS = 24
NB_CPUS = 2


@pytest.fixture(scope="module", autouse=True)
def tempdir():
    tempdir_path = Path("utility_files/tmp").resolve()
    if tempdir_path.is_dir():
        shutil.rmtree(tempdir_path)
    tempdir_path.mkdir()
    yield tempdir_path
    # all code after yield is executed as teardown
    print("Teardown: removing tmp directory")
    shutil.rmtree(tempdir_path)


# =============== GRAVITO ACOUSTIC TEST ===============
@pytest.fixture(scope="module")
def dfs_gravito_acoustic(tempdir):
    config = {
        "equilibrium_type": "gravito_acoustic",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": np.linspace(0, np.sqrt(250), NB_RUNS),
            "k3": np.linspace(0, np.sqrt(250), NB_RUNS),
            "cte_p0": 1,
            "g": 0.5,
            "alpha": 20.42,
        },
        "basename_datfile": "gravito_acoustic",
        "output_folder": str(tempdir),
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    }
    parfiles = pylbo.generate_parfiles(config, output_dir=tempdir)
    pylbo.run_legolas(parfiles, remove_parfiles=True, nb_cpus=NB_CPUS)
    datfiles = sorted(Path(tempdir).glob("*gravito_acoustic.dat"))
    if not datfiles:
        raise ValueError("No datfiles found!")
    return datfiles


@pytest.mark.mpl_image_compare(
    baseline_dir=BASELINE_DIR,
    filename="gravito_acoustic.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_gravito_acoustic(dfs_gravito_acoustic):
    datasets = pylbo.load(dfs_gravito_acoustic)
    ms = pylbo.MultiSpectrum(datasets)
    fig, ax = ms.plot_precoded_run()
    return fig


# =============== GRAVITO-MHD TEST ===============
@pytest.fixture(scope="module")
def dfs_gravito_mhd(tempdir):
    config = {
        "equilibrium_type": "gravito_mhd",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": np.linspace(0, np.sqrt(250), NB_RUNS),
            "k3": np.linspace(0, np.sqrt(250), NB_RUNS),
            "cte_p0": 0.5,
            "g": 20,
            "alpha": 20,
        },
        "basename_datfile": "gravito_mhd",
        "output_folder": str(tempdir),
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    }
    parfiles = pylbo.generate_parfiles(config, output_dir=tempdir)
    pylbo.run_legolas(parfiles, remove_parfiles=True, nb_cpus=NB_CPUS)
    datfiles = sorted(Path(tempdir).glob("*gravito_mhd.dat"))
    if not datfiles:
        raise ValueError("No datfiles found!")
    return datfiles


@pytest.mark.mpl_image_compare(
    baseline_dir=BASELINE_DIR,
    filename="gravito_mhd.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_gravito_mhd(dfs_gravito_mhd):
    datasets = pylbo.load(dfs_gravito_mhd)
    ms = pylbo.MultiSpectrum(datasets)
    fig, ax = ms.plot_precoded_run()
    return fig


# =============== INTERCHANGE MODES TEST ===============
@pytest.fixture(scope="module")
def dfs_interchange(tempdir):
    config = {
        "equilibrium_type": "interchange_modes",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": np.pi * np.sin(np.linspace(0, np.pi, NB_RUNS)),
            "k3": np.pi * np.cos(np.linspace(0, np.pi, NB_RUNS)),
            "g": 0.5,
            "cte_p0": 0.25,
            "lambda": 0.3,
            "alpha": 20.0,
        },
        "basename_datfile": "interchange",
        "output_folder": str(tempdir),
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    }
    parfiles = pylbo.generate_parfiles(config, output_dir=tempdir)
    pylbo.run_legolas(parfiles, remove_parfiles=True, nb_cpus=NB_CPUS)
    datfiles = sorted(Path(tempdir).glob("*interchange.dat"))
    if not datfiles:
        raise ValueError("No datfiles found!")
    return datfiles


@pytest.mark.mpl_image_compare(
    baseline_dir=BASELINE_DIR,
    filename="interchange.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_interchange(dfs_interchange):
    datasets = pylbo.load(dfs_interchange)
    ms = pylbo.MultiSpectrum(datasets)
    fig, ax = ms.plot_precoded_run(annotate_continua=True)
    return fig


# =============== CONSTANT CURRENT TEST ===============
@pytest.fixture(scope="module")
def dfs_constant_current(tempdir):
    config = {
        "equilibrium_type": "constant_current_tokamak",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": -2.0,
            "k3": 0.2,
            "j0": (2.0 * 0.2) / np.linspace(1.9, 2.1, NB_RUNS),
            "cte_rho0": 1.0,
            "cte_B03": 1.0,
        },
        "basename_datfile": "constant_current",
        "output_folder": str(tempdir),
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    }
    parfiles = pylbo.generate_parfiles(config, output_dir=tempdir)
    pylbo.run_legolas(parfiles, remove_parfiles=False, nb_cpus=NB_CPUS)
    datfiles = sorted(Path(tempdir).glob("*constant_current.dat"))
    if not datfiles:
        raise ValueError("No datfiles found!")
    return datfiles


@pytest.mark.mpl_image_compare(
    baseline_dir=BASELINE_DIR,
    filename="constant_current.png",
    savefig_kwargs={"dpi": FIG_DPI},
    tolerance=10,
)
def test_constant_current(dfs_constant_current):
    datasets = pylbo.load(dfs_constant_current)
    ms = pylbo.MultiSpectrum(datasets)
    fig, ax = ms.plot_precoded_run(annotate_continua=True)
    return fig


# =============== PHOTOSPHERIC FLUXTUBE TEST ===============
@pytest.fixture(scope="module")
def dfs_photospheric(tempdir):
    config = {
        "equilibrium_type": "photospheric_flux_tube",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": 0.0,
            "k3": np.linspace(0.1, 6.2, NB_RUNS),
            "cte_rho0": 1.0,
            "cte_p0": 1.0,
            "r0": 1.0,
        },
        "mesh_accumulation": True,
        "basename_datfile": "photospheric_tube",
        "output_folder": str(tempdir),
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    }
    parfiles = pylbo.generate_parfiles(config, output_dir=tempdir)
    pylbo.run_legolas(parfiles, remove_parfiles=True, nb_cpus=NB_CPUS)
    datfiles = sorted(Path(tempdir).glob("*photospheric_tube.dat"))
    if not datfiles:
        raise ValueError("No datfiles found!")
    return datfiles


@pytest.mark.mpl_image_compare(
    baseline_dir=BASELINE_DIR,
    filename="photospheric_tube.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_photospheric(dfs_photospheric):
    datasets = pylbo.load(dfs_photospheric)
    ms = pylbo.MultiSpectrum(datasets)
    fig, ax = ms.plot_precoded_run()
    return fig


# =============== CORONAL FLUXTUBE TEST ===============
@pytest.fixture(scope="module")
def dfs_coronal(tempdir):
    config = {
        "equilibrium_type": "coronal_flux_tube",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": 0.0,
            "k3": np.linspace(0.1, 6.2, NB_RUNS),
            "cte_rho0": 1.0,
            "cte_p0": 1.0,
            "r0": 1.0,
        },
        "mesh_accumulation": True,
        "basename_datfile": "coronal_tube",
        "output_folder": str(tempdir),
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    }
    parfiles = pylbo.generate_parfiles(config, output_dir=tempdir)
    pylbo.run_legolas(parfiles, remove_parfiles=True, nb_cpus=NB_CPUS)
    datfiles = sorted(Path(tempdir).glob("*coronal_tube.dat"))
    if not datfiles:
        raise ValueError("No datfiles found!")
    return datfiles


@pytest.mark.mpl_image_compare(
    baseline_dir=BASELINE_DIR,
    filename="coronal_tube.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_coronal(dfs_coronal):
    datasets = pylbo.load(dfs_coronal)
    ms = pylbo.MultiSpectrum(datasets)
    fig, ax = ms.plot_precoded_run()
    return fig
