import numpy as np
import pytest
from pathlib import Path
import pylbo

# =============== ADIABATIC TESTS ===============


@pytest.fixture(scope="module", autouse=True)
def ds_adiab():
    file = Path("utility_files/v1_adiabatic_car.dat").resolve()
    ds = pylbo.load(file)
    return ds


@pytest.fixture(scope="module", autouse=True)
def ds_adiab_cyl():
    file = Path("utility_files/v1_adiabatic_cyl.dat").resolve()
    ds = pylbo.load(file)
    return ds


def test_cs_adiab(ds_adiab):
    cs = ds_adiab.get_sound_speed()
    assert np.all(cs == pytest.approx(np.sqrt(5 / 3)))


def test_cs_adiab_equal(ds_adiab):
    cs = ds_adiab.get_sound_speed()
    # should be equal everywhere
    assert np.all(np.diff(cs) == pytest.approx(0))


def test_ca_adiab(ds_adiab):
    ca = ds_adiab.get_alfven_speed()
    assert np.all(ca == pytest.approx(1))


def test_ca_adiab_equal(ds_adiab):
    ca = ds_adiab.get_alfven_speed()
    assert np.all(np.diff(ca) == pytest.approx(0))


def test_ct_adiab_cart(ds_adiab):
    pylbo.set_loglevel("error")
    ct = ds_adiab.get_tube_speed()
    assert ct is None
    pylbo.set_loglevel("warning")


def test_ct_adiab_cyl(ds_adiab_cyl):
    ct = ds_adiab_cyl.get_tube_speed()
    assert np.all(ct == pytest.approx(np.sqrt(5 / 3) / np.sqrt(8 / 3)))


def test_ct_adiab_cyl_equal(ds_adiab_cyl):
    ct = ds_adiab_cyl.get_tube_speed()
    assert np.all(np.diff(ct) == pytest.approx(0))


def test_k0sq_adiab(ds_adiab):
    k0sq = ds_adiab.get_k0_squared()
    assert k0sq == pytest.approx(np.pi ** 2)


def test_re_adiab(ds_adiab):
    pylbo.set_loglevel("error")
    re = ds_adiab.get_reynolds()
    assert re is None
    pylbo.set_loglevel("warning")


def test_rm_adiab(ds_adiab):
    pylbo.set_loglevel("error")
    rm = ds_adiab.get_reynolds_magnetic()
    assert rm is None
    pylbo.set_loglevel("warning")


def test_nearest_evs():
    file = Path("utility_files/v1_datfile_efs.dat").resolve()
    ds = pylbo.load(file)
    guesses = [11.1851, -2.67131, -2.48601]
    idxs, evs = ds.get_nearest_eigenvalues(guesses)
    assert evs == pytest.approx(np.array([11.1851013, -2.67130721, -2.48601357]))
    assert idxs == pytest.approx([36, 43, 64])


# =============== RESISTIVE TESTS ===============


@pytest.fixture(scope="module", autouse=True)
def ds_resis():
    file = Path("utility_files/v1_resistive.dat").resolve()
    ds = pylbo.load(file)
    return ds


def test_cs_resis_equal(ds_resis):
    cs = ds_resis.get_sound_speed()
    assert np.all(np.diff(cs) == pytest.approx(0))


def test_cs_resis(ds_resis):
    cs = ds_resis.get_sound_speed()
    assert np.all(cs == pytest.approx(np.sqrt(0.375 / 3)))


def test_ca_resis_equal(ds_resis):
    ca = ds_resis.get_alfven_speed()
    assert np.all(ca == pytest.approx(1))


def test_ca_resis(ds_resis):
    ca = ds_resis.get_alfven_speed()
    assert np.all(ca == pytest.approx(1))


def test_k0sq_resis(ds_resis):
    k0sq = ds_resis.get_k0_squared()
    assert k0sq == pytest.approx(0.49 ** 2)


def test_re_resis_equal(ds_resis):
    re = ds_resis.get_reynolds()
    assert np.all(np.diff(re) == pytest.approx(0))


def test_re_resis(ds_resis):
    re = ds_resis.get_reynolds()
    assert np.all(re == pytest.approx(np.sqrt(0.375 / 3) * 1e4))


def test_rm_resis_equal(ds_resis):
    rm = ds_resis.get_reynolds_magnetic()
    assert np.all(np.diff(rm) == pytest.approx(0, abs=1e-10))


def test_rm_resis(ds_resis):
    rm = ds_resis.get_reynolds_magnetic()
    assert np.all(rm == pytest.approx(1e4))
