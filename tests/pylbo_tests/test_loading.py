import pytest
import pylbo


def test_invalid_file():
    with pytest.raises(FileNotFoundError):
        pylbo.load("unknown_file")


def test_invalid_suffix(tmpdir):
    from pylbo.exceptions import InvalidLegolasFile

    parfile = tmpdir / "parfile.par"
    parfile.write_text("content")
    with pytest.raises(InvalidLegolasFile):
        pylbo.load(parfile)


def test_load_invalid(datv0):
    with pytest.raises(ValueError):
        pylbo.load([datv0, datv0])


def test_load_series_empty():
    with pytest.raises(ValueError):
        pylbo.load_series([])


@pytest.mark.timeout(5)
def test_load_logfile_v0(logv0):
    pylbo.load_logfile(logv0)


@pytest.mark.timeout(5)
def test_load_datfile_v0(datv0):
    pylbo.load(datv0)


@pytest.mark.timeout(5)
def test_load_datfile_v090(datv090):
    pylbo.load(datv090)


@pytest.mark.timeout(5)
def test_load_datfile_v100(datv100):
    pylbo.load(datv100)


@pytest.mark.timeout(5)
def test_load_datfile_v112(datv112):
    pylbo.load(datv112)


@pytest.mark.timeout(5)
def test_load_series(datv100):
    series = pylbo.load_series([datv100, datv100, datv100])
    assert len(series) == 3
