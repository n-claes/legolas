import pytest
import shutil
from pathlib import Path
import pylbo

pylbo.set_loglevel("warning")

imagedirpath = Path(__file__).parent / "image_comparisons"
KEEP_FILES_OPTION = "--keep-files"

# @note: all fixtures defined here will be accessible to all tests
#        in the same directory as this file.


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    def remove_results_dir():
        resultsdirpath = (Path(__file__).parent / "results").resolve()
        if resultsdirpath.is_dir():
            if not request.config.getoption(KEEP_FILES_OPTION):
                shutil.rmtree(resultsdirpath)

    def remove_image_dir():
        # only remove this directory if contents are empty (all tests passed)
        try:
            imagedirpath.rmdir()
        except OSError:
            pass

    request.addfinalizer(remove_results_dir)
    request.addfinalizer(remove_image_dir)


def pytest_addoption(parser):
    parser.addoption(
        KEEP_FILES_OPTION,
        action="store_true",
        help="if supplied, does not remove files after test completion.",
    )


@pytest.fixture()
def keep_files(request):
    return request.config.getoption(KEEP_FILES_OPTION)


@pytest.fixture
def imagedir(scope="session", autouse=True):
    if not imagedirpath.is_dir():
        imagedirpath.mkdir()
    yield imagedirpath


@pytest.fixture
def ds_test(setup):
    if setup["test_needs_run"]:
        parfile = pylbo.generate_parfiles(
            parfile_dict=setup["config"],
            basename=setup["datfile"].stem,
            output_dir=setup["config"]["output_folder"],
        )
        pylbo.run_legolas(parfile, remove_parfiles=True)
        setup["test_needs_run"] = False
    return pylbo.load(setup["datfile"])


@pytest.fixture
def log_test(setup):
    return pylbo.load_logfile(setup["logfile"], sort=True)


@pytest.fixture
def ds_answer(setup):
    return pylbo.load(setup["answer_datfile"])


@pytest.fixture
def log_answer(setup):
    return pylbo.load_logfile(setup["answer_logfile"], sort=True)


# ===== FIXTURES FOR MULTIRUNS =====
@pytest.fixture()
def tempdir():
    tempdir_path = (Path(__file__).parent / "results/tmp").resolve()
    if tempdir_path.is_dir():
        shutil.rmtree(tempdir_path)
    tempdir_path.mkdir()
    yield tempdir_path
    # all code after yield is executed as teardown
    shutil.rmtree(tempdir_path)


@pytest.fixture
def series_test(tempdir, setup):
    nb_cpus = 2
    setup["config"]["output_folder"] = str(tempdir)
    parfiles = pylbo.generate_parfiles(
        parfile_dict=setup["config"],
        basename=setup["name"],
        output_dir=tempdir,
    )
    print("")
    pylbo.run_legolas(parfiles, remove_parfiles=True, nb_cpus=nb_cpus)
    datfiles = sorted(tempdir.glob(f"*{setup['config']['basename_datfile']}.dat"))
    return pylbo.load_series(datfiles)
