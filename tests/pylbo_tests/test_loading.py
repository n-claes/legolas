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
