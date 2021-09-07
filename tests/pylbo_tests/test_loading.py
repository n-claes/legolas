import pytest
import numpy as np
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


def test_load_no_version(datv0):
    ds = pylbo.load(datv0)
    assert ds.legolas_version == "0.0.0"


def test_load_series_empty():
    with pytest.raises(ValueError):
        pylbo.load_series([])


def test_load_multiple_equilibria(datv1, datv112_eta):
    with pytest.raises(ValueError):
        pylbo.load_series([datv1, datv112_eta])


def test_load_logfile(logv0):
    eigenvals = pylbo.load_logfile(logv0)
    assert isinstance(eigenvals, np.ndarray)


def test_load_logfile_and_sort(logv0):
    eigenvals = pylbo.load_logfile(logv0, sort=True)
    assert isinstance(eigenvals, np.ndarray)
    assert np.all(eigenvals[:-1] <= eigenvals[1:])
