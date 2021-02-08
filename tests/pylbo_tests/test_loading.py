import pytest
from pathlib import Path
import pylbo

utils = (Path(__file__).parent / "utility_files").resolve()


@pytest.mark.timeout(10)
@pytest.fixture(scope="module")
def logfile_v0():
    return pylbo.load_logfile(utils / "v0_logfile_efs.log")


@pytest.mark.timeout(10)
@pytest.fixture(scope="module")
def datfile_v0():
    return pylbo.load(utils / "v0_datfile_efs.dat")


@pytest.mark.timeout(10)
@pytest.fixture(scope="module")
def datfile_v09():
    return pylbo.load(utils / "v0.9.0_datfile.dat")


@pytest.mark.timeout(10)
@pytest.fixture(scope="module")
def datfile_v1():
    return pylbo.load(utils / "v1_datfile_matrices.dat")


@pytest.mark.timeout(10)
@pytest.fixture(scope="module")
def series(datfile_v1):
    return pylbo.load_series([datfile_v1, datfile_v1, datfile_v1])


@pytest.mark.timeout(10)
def test_invalid_file():
    with pytest.raises(FileNotFoundError):
        pylbo.load("unknown_file")


def test_invalid_suffix(tmp_path):
    from pylbo.exceptions import InvalidLegolasFile

    parfile = tmp_path / "parfile.par"
    parfile.write_text("content")
    with pytest.raises(InvalidLegolasFile):
        pylbo.load(parfile)


def test_load_invalid(datfile_v0, datfile_v1):
    with pytest.raises(ValueError):
        pylbo.load([datfile_v0, datfile_v1])


def test_load_series_empty():
    with pytest.raises(ValueError):
        pylbo.load_series([])
