import pytest
import numpy as np
from pathlib import Path
import pylbo

pylbo.set_loglevel("warning")


@pytest.mark.timeout(10)
def test_load_v0():
    file = Path("utility_files/v0_datfile_efs.dat").resolve()
    pylbo.load(file)


@pytest.mark.timeout(10)
def test_load_v1():
    file = Path("utility_files/v1_datfile_efs.dat").resolve()
    pylbo.load(file)


def test_filenotfound():
    file = "non-existent-file"
    with pytest.raises(FileNotFoundError):
        pylbo.load(file)


def test_filenotdat():
    from pylbo.utilities.exceptions import InvalidLegolasFile

    file = Path("utility_files/v0_logfile_efs.log").resolve()
    with pytest.raises(InvalidLegolasFile):
        pylbo.load(file)


def test_load_multi():
    from pylbo.data_management.data_container import LegolasDataContainer

    files = [
        Path("utility_files/v0_datfile_efs.dat").resolve(),
        Path("utility_files/v1_datfile_efs.dat").resolve(),
    ]
    datasets = pylbo.load(files)
    for ds in datasets:
        assert isinstance(ds, LegolasDataContainer)


def test_load_listoneitem():
    from pylbo.data_management.data_container import LegolasDataContainer

    file = [Path("utility_files/v1_datfile_efs.dat").resolve()]
    ds = pylbo.load(file)
    assert isinstance(ds, LegolasDataContainer)


def test_load_invalidlog():
    file = Path("utility_files/v1_datfile_efs.dat").resolve()
    with pytest.raises(ValueError):
        pylbo.read_log_file(file)


def test_load_nonexistentlog():
    file = "non-existent-file.log"
    with pytest.raises(FileNotFoundError):
        pylbo.read_log_file(file)


def test_load_log():
    file = Path("utility_files/v0_logfile_efs.log").resolve()
    evs = pylbo.read_log_file(file)
    assert isinstance(evs, np.ndarray)
    for ev in evs:
        assert isinstance(ev, complex)


def test_load_log_sort():
    file = Path("utility_files/v0_logfile_efs.log").resolve()
    evs = pylbo.read_log_file(file, sort=True)
    assert isinstance(evs, np.ndarray)
    for ev in evs:
        assert isinstance(ev, complex)
    assert np.all(np.diff(evs) >= 0)
