import pylbo.gimli as gimli
import sympy as sp
import numpy as np
from scipy.io import FortranFile
import filecmp


def test_variables():
    obj = gimli.Variables()
    keychain = obj.__dict__.keys()
    for key in [
        "x",
        "y",
        "z",
        "rho0",
        "T0",
        "B0sq",
        "k2",
        "k3",
        "rhoc",
        "Tc",
        "B2c",
        "B3c",
        "v2c",
        "v3c",
        "pc",
        "p1",
        "p2",
        "p3",
        "p4",
        "p5",
        "p6",
        "p7",
        "p8",
        "alpha",
        "beta",
        "delta",
        "theta",
        "tau",
        "lam",
        "nu",
        "r0",
        "rc",
        "rj",
        "Bth0",
        "Bz0",
        "V",
        "j0",
        "g",
    ]:
        assert key in keychain
        assert isinstance(obj.__dict__[key], sp.Symbol)
    assert "fkey" in keychain
    assert isinstance(obj.__dict__["fkey"], dict)


def test_equilibrium():
    var = gimli.Variables()
    obj = gimli.Equilibrium(var, var.rhoc, 0, 0, var.Tc)
    keychain = obj.__dict__.keys()
    for key in ["variables", "rho0", "v02", "v03", "T0", "B02", "B03", "_dict_phys"]:
        assert key in keychain
    assert isinstance(obj.get_physics(), dict)
    assert isinstance(obj.get_dependencies(), dict)


def test_legolas_userfile_hd(tmpdir, mod_usr_hd):
    config = {
        "geometry": "Cartesian",
        "x_start": -1,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {
            "k2": 1,
            "k3": 0,
            "cte_rho0": 1,
            "cte_T0": 0.5,
        },
        "equilibrium_type": "user_defined",
        "physics_type": "hd",
        "logging_level": 1,
    }
    var = gimli.Variables()
    eq = gimli.Equilibrium(var, var.rhoc, 0, 0, var.Tc)
    obj = gimli.Legolas(eq, config)
    obj.user_module(filename="smod_user_hd", loc=tmpdir)
    assert filecmp.cmp(
        str((tmpdir / "smod_user_hd.f08").resolve()), str(mod_usr_hd), shallow=False
    )


def test_legolas_userfile_mhd(tmpdir, mod_usr_mhd):
    config = {
        "geometry": "Cartesian",
        "x_start": -1,
        "x_end": 1,
        "gridpoints": 51,
        "parameters": {
            "k2": 1,
            "k3": 0,
            "cte_rho0": 1,
            "cte_T0": 0.5,
            "cte_B02": 0.25,
        },
        "equilibrium_type": "user_defined",
        "physics_type": "mhd",
        "logging_level": 1,
    }
    var = gimli.Variables()
    eq = gimli.Equilibrium(var, var.rhoc, 0, 0, var.Tc, B02=var.B2c, B03=0)
    obj = gimli.Legolas(eq, config)
    obj.user_module(filename="smod_user_mhd", loc=tmpdir)
    assert filecmp.cmp(
        str((tmpdir / "smod_user_mhd.f08").resolve()), str(mod_usr_mhd), shallow=False
    )


def test_amrvac_preparation(tmpdir, datv211_harris, vacv211_harris):
    config = {
        "datfile": datv211_harris,
        "physics_type": "mhd",
        "ev_guess": [0.01636j, 1.397e-2 - 2.843e-4 * 1j, -1.397e-2 - 2.843e-4 * 1j],
        "percentage": 0.01,
        "quantity": "B02",
    }
    amrvac = gimli.Amrvac(config)
    amrvac.prepare_legolas_data(loc=tmpdir)

    base = FortranFile(vacv211_harris, "r")
    test = FortranFile(tmpdir / "v2.1.1_harris.ldat", "r")
    for ii in range(2):
        base_data = base.read_ints(dtype=np.int32)
        test_data = test.read_ints(dtype=np.int32)
        assert np.array_equal(base_data, test_data)

    for ii in range(11):
        base_data = base.read_reals(dtype=np.float64)
        test_data = test.read_reals(dtype=np.float64)
        assert np.allclose(base_data, test_data, rtol=1e-8, atol=1e-10)


def test_numerical_equilibrium(tmpdir, numerical_lar):
    x = np.linspace(-np.pi, np.pi, 1000)
    dictionary = {"x": x, "rho0": 2.0 + np.sin(x), "T0": 1.0 / (2.0 + np.sin(x))}
    equil = gimli.NumericalEquilibrium(dictionary)
    equil.to_legolas_arrays(filename="test_numerical", loc=str(tmpdir.resolve()))

    base = FortranFile(numerical_lar, "r")
    test = FortranFile(tmpdir / "test_numerical.lar", "r")
    base_data = base.read_ints(dtype=np.int32)
    test_data = test.read_ints(dtype=np.int32)
    assert np.array_equal(base_data, test_data)

    for ii in range(10):
        base_data = base.read_reals(dtype=np.float64)
        test_data = test.read_reals(dtype=np.float64)
        assert np.allclose(base_data, test_data, rtol=1e-8, atol=1e-10)
