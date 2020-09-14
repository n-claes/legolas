import pytest
import copy
import numpy as np
from pathlib import Path
import pylbo
from pylbo.utilities.continua import _thermal_continuum

pylbo.set_loglevel("warning")


@pytest.fixture(scope="module", autouse=True)
def ds_adiab():
    file = Path("utility_files/v1_adiabatic_car.dat").resolve()
    return pylbo.load(file)


@pytest.fixture(scope="module", autouse=True)
def ds_thermal():
    file = Path("utility_files/v1_datfile_magnetothermal.dat").resolve()
    return pylbo.load(file)


def test_adiab(ds_adiab):
    wth = _thermal_continuum(ds_adiab, np.ones_like(ds_adiab.grid_gauss))
    assert np.all(wth == pytest.approx(0))


def test_tempzero(ds_thermal):
    copyds = copy.deepcopy(ds_thermal)
    T0 = copyds.equilibria.get("T0")
    copyds.equilibria.update({"T0": np.zeros_like(T0)})
    wth = _thermal_continuum(copyds, copyds.grid_gauss)
    assert np.all(wth == pytest.approx(0))


def test_wthermal(ds_thermal):
    wth = ds_thermal.continua.get("wth")
    assert np.min(wth) == pytest.approx(0.002793501)
    assert np.max(wth) == pytest.approx(0.013565217)


def test_wslow(ds_thermal):
    wspos = ds_thermal.continua.get("wS+")
    wsneg = ds_thermal.continua.get("wS-")
    assert np.all(wspos == pytest.approx(0))
    assert np.all(wsneg == pytest.approx(0))


def test_walfven(ds_thermal):
    wapos = ds_thermal.continua.get("wA+")
    waneg = ds_thermal.continua.get("wA-")
    assert np.all(wapos == pytest.approx(0))
    assert np.all(waneg == pytest.approx(0))


def test_doppler(ds_thermal):
    doppler = ds_thermal.continua.get("dopp")
    assert np.all(doppler == pytest.approx(0))



