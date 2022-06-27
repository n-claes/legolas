import os
from pathlib import Path
import shutil
import pytest
import pylbo

pylbo.set_loglevel("warning")

KEEP_FILES_OPTION = "--keep-files"
GENERATE_BASELINE_OPTION = "--generate"


def check_folder(folder):
    if not folder.is_dir():
        folder.mkdir()


def safely_remove_tree(root_folder):
    # walks through the root folder, removes directories only if they are empty
    folders = list(os.walk(root_folder))
    for folder in reversed(folders):
        try:
            Path(folder[0]).rmdir()
        except OSError:
            pass


@pytest.fixture(scope="session", autouse=True)
def maindir():
    return Path(__file__).parent.resolve()


@pytest.fixture(scope="session", autouse=True)
def baselinedir(maindir):
    return maindir / "baseline"


@pytest.fixture(scope="session", autouse=True)
def resultsdir(maindir):
    folder = maindir / "test_results"
    check_folder(folder)
    return folder


@pytest.fixture(scope="class", autouse=True)
def casedir(request, resultsdir):
    folder = resultsdir / request.cls.filename
    check_folder(folder)
    return folder


@pytest.fixture(scope="class", autouse=True)
def datfiledir(casedir):
    folder = casedir / "datfiles"
    check_folder(folder)
    return folder


@pytest.fixture(scope="class", autouse=True)
def spectradir(casedir):
    folder = casedir / "spectra"
    check_folder(folder)
    return folder


@pytest.fixture(scope="class", autouse=True)
def eigfuncdir(casedir):
    folder = casedir / "eigenfunctions"
    check_folder(folder)
    return folder


@pytest.fixture(scope="session", autouse=True)
def keep_files(request):
    return request.config.getoption(KEEP_FILES_OPTION)


@pytest.fixture(scope="class", autouse=True)
def set_directories_as_cls_attrs(
    request, keep_files, resultsdir, casedir, datfiledir, spectradir, eigfuncdir
):
    request.cls._keep_files = keep_files
    request.cls._resultsdir = resultsdir
    request.cls._datfiledir = datfiledir
    request.cls._casedir = casedir
    request.cls._spectradir = spectradir
    request.cls._eigfuncdir = eigfuncdir


def pytest_addoption(parser):
    parser.addoption(
        KEEP_FILES_OPTION,
        action="store_true",
        help="if supplied, does not remove files after test completion.",
    )
    parser.addoption(
        GENERATE_BASELINE_OPTION,
        action="store_true",
        dest="generate_baseline",
        help="if set to true, (re)generates the baseline dataset.",
    )


@pytest.fixture(scope="class", autouse=True)
def case_cleanup(request, casedir, datfiledir):
    def remove_datfiledir():
        # always remove generated datfiles by default
        if not request.config.getoption(KEEP_FILES_OPTION):
            shutil.rmtree(datfiledir)

    def remove_casedir():
        safely_remove_tree(casedir)

    request.addfinalizer(remove_datfiledir)
    request.addfinalizer(remove_casedir)


@pytest.fixture(scope="session", autouse=True)
def final_cleanup(request, resultsdir):
    def remove_resultsdir():
        safely_remove_tree(resultsdir)

    request.addfinalizer(remove_resultsdir)


def pytest_make_parametrize_id(val):
    if isinstance(val, dict):
        # spectrum images
        if "xlim" in val.keys():
            return f"Re(w)={val['xlim']}, Im(w)={val['ylim']}"
        # eigenfunction images
        if "eigenvalue" in val.keys():
            return f"w={val['eigenvalue']:.8f}"
    return None
