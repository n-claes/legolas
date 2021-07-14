import pytest
import shutil
from pathlib import Path
import pylbo
import pylbo.testing

pylbo.set_loglevel("warning")

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
def tmpdir(scope="session", autouse=True):
    if not tmpdir_path.is_dir():
        tmpdir_path.mkdir()
    yield tmpdir_path


@pytest.fixture
def logv0():
    return utils / "v0_logfile_efs.log"


@pytest.fixture
def datv0():
    return utils / "v0_datfile_efs.dat"


@pytest.fixture
def datv090():
    return utils / "v0.9.0_datfile.dat"


@pytest.fixture
def datv100():
    return utils / "v1_datfile_matrices.dat"


@pytest.fixture
def datv112():
    return utils / "v1.1.2_datfile_efs.dat"


@pytest.fixture
def ds_v100(datv100):
    return pylbo.load(datv100)


@pytest.fixture
def ds_v112(datv112):
    return pylbo.load(datv112)
