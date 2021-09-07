import pytest
import shutil
import matplotlib.pyplot as plt
from pathlib import Path
import pylbo
import pylbo.testing
from pylbo.visualisation.continua import ContinuaHandler


pylbo.set_loglevel("error")

KEEP_FILES_OPTION = "--keep-files"
tmpdir_path = Path(__file__).resolve().parent / "tmp"
utils = Path(__file__).resolve().parent / "utility_files"


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    def remove_tmp_dir():
        if tmpdir_path.is_dir():
            if not request.config.getoption(KEEP_FILES_OPTION):
                shutil.rmtree(tmpdir_path)

    request.addfinalizer(remove_tmp_dir)


@pytest.fixture(autouse=True)
def close_figures_after_test():
    yield
    plt.close("all")


def pytest_addoption(parser):
    parser.addoption(
        KEEP_FILES_OPTION,
        action="store_true",
        help="if supplied, does not remove files after test completion.",
    )


@pytest.fixture
def keep_files(request):
    return request.config.getoption(KEEP_FILES_OPTION)


@pytest.fixture
def tmpdir(scope="session"):
    if not tmpdir_path.is_dir():
        tmpdir_path.mkdir()
    yield tmpdir_path


@pytest.fixture
def default_pf_dict(tmpdir):
    config = {
        "gridpoints": 10,
        "equilibrium_type": "suydam_cluster",
        "show_results": False,
        "write_eigenfunctions": False,
        "basename_datfile": "default_ds",
        "output_folder": str(tmpdir),
        "logging_level": 0,
    }
    return config


@pytest.fixture
def default_parfile(tmpdir, default_pf_dict):
    return pylbo.generate_parfiles(default_pf_dict, output_dir=tmpdir)


@pytest.fixture
def default_ds(tmpdir, default_pf_dict):
    filepath = (tmpdir / "default_ds.dat").resolve()
    if not filepath.is_file():
        parfile = pylbo.generate_parfiles(default_pf_dict, output_dir=tmpdir)
        pylbo.run_legolas(parfile)
    return pylbo.load(filepath)


@pytest.fixture
def fake_ds():
    datfile = utils / "v1.1.2_datfile_efs.dat"
    return pylbo.testing.FakeDataSet(datfile=datfile, seed=20210715)


@pytest.fixture
def c_handle():
    return ContinuaHandler(interactive=True)


@pytest.fixture
def logv0():
    return utils / "v0_logfile_efs.log"


@pytest.fixture
def datv0():
    return utils / "v0_datfile_efs.dat"


@pytest.fixture
def datv1():
    return utils / "v1_datfile_matrices.dat"


@pytest.fixture
def datv112_eta():
    return utils / "v1.1.2_datfile_eta.dat"


@pytest.mark.timeout(5)
@pytest.fixture
def ds_v090():
    return pylbo.load(utils / "v0.9.0_datfile.dat")


@pytest.mark.timeout(5)
@pytest.fixture
def ds_v100():
    return pylbo.load(utils / "v1_datfile_matrices.dat")


@pytest.fixture
def ds_v112():
    return pylbo.load(utils / "v1.1.2_datfile_efs.dat")


@pytest.mark.timeout(5)
@pytest.fixture
def ds_v112_eta():
    return pylbo.load(utils / "v1.1.2_datfile_eta.dat")


@pytest.mark.timeout(5)
@pytest.fixture
def ds_v114_subset():
    return pylbo.load(utils / "v1.1.4_datfile_subset.dat")


@pytest.mark.timeout(5)
@pytest.fixture
def ds_v114_subset_defs():
    return pylbo.load(utils / "v1.1.4_datfile_subset_defs.dat")


@pytest.mark.timeout(5)
@pytest.fixture
def ds_v114():
    return pylbo.load(utils / "v1.1.4_datfile.dat")


@pytest.mark.timeout(5)
@pytest.fixture
def series_v100():
    return pylbo.load_series([utils / "v1_datfile_matrices.dat"] * 3)


@pytest.mark.timeout(5)
@pytest.fixture
def series_v112():
    return pylbo.load_series([utils / "v1.1.2_datfile_efs.dat"] * 3)


@pytest.mark.timeout(5)
@pytest.fixture
def series_v112_eta():
    return pylbo.load_series([utils / "v1.1.2_datfile_eta.dat"] * 5)
