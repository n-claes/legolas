import pytest
import shutil
from pathlib import Path
import pylbo

pylbo.set_loglevel("warning")

# @note: all fixtures defined here will be accessible to all tests
#        in the same directory as this file.


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
def eigfuncs_test(ds_test, setup):
    if setup.get("ev_guesses", None) is not None:
        return ds_test.get_eigenfunctions(ev_guesses=setup["ev_guesses"])
    else:
        return None


@pytest.fixture
def ds_answer(setup):
    return pylbo.load(setup["answer_datfile"])


@pytest.fixture
def log_answer(setup):
    return pylbo.load_logfile(setup["answer_logfile"], sort=True)


@pytest.fixture
def eigfuncs_answer(ds_answer, setup):
    if setup.get("ev_guesses", None) is not None:
        return ds_answer.get_eigenfunctions(ev_guesses=setup["ev_guesses"])
    else:
        return None


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
