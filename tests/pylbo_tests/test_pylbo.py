import pytest
import pylbo
import numpy as np
from pathlib import Path

output = (pylbo.LEGOLAS_DIR / 'tests/pylbo_tests').resolve()
datfile = (output / 'resistive_homo.dat').resolve()
logfile = (output / 'resistive_homo.log').resolve()

config = {
    'gridpoints': 51,
    'geometry': 'Cartesian',
    'x_start': 0,
    'x_end': 1,
    'equilibrium_type': 'resistive_homo',
    'logging_level': 0,
    'show_results': False,
    'write_eigenfunctions': True,
    'write_matrices': True,
    'basename_datfile': datfile.stem,
    'basename_logfile': logfile.stem,
    'output_folder': str(output)
}

@pytest.fixture(scope='module', autouse=True)
def ds():
    parfile = pylbo.generate_parfiles(parfile_dict=config, basename_parfile=datfile.stem,
                                      output_dir=output)
    pylbo.run_legolas(parfile, remove_parfiles=True)
    return pylbo.load(datfile)

@pytest.fixture(scope='module', autouse=True)
def log():
    return pylbo.read_log_file(logfile, sort=True)

@pytest.fixture(scope='module', autouse=True)
def header():
    with open(datfile, 'rb') as istream:
        hdr = pylbo.get_header(istream)
    return hdr

def test_existence():
    assert Path.is_file(datfile)
    assert Path.is_file(logfile)
    assert Path.is_dir(output)

def test_geometry(ds):
    assert ds.geometry == 'Cartesian'

def test_x_start(ds):
    assert ds.x_start == 0

def test_x_end(ds):
    assert ds.x_end == 1

def test_gridpts(ds):
    gridpts = 51
    assert ds.gridpts == gridpts
    assert ds.gauss_gridpts == 4 * (gridpts - 1)
    assert ds.matrix_gridpts == 16 * gridpts
    assert ds.ef_gridpts == 2 * gridpts - 1

def test_gamma(ds):
    assert ds.gamma == 5/3

def test_eq_type(ds):
    assert ds.eq_type == 'resistive_homo'

def test_grid(ds):
    assert ds.grid[0] == 0
    assert ds.grid[-1] == 1

def test_evs(ds, log):
    with open(datfile, 'rb') as istream:
        header = pylbo.get_header(istream)
        evs_dat = pylbo.read_eigenvalues(istream, header, omit_large_evs=False)

    evs_dat = np.sort(evs_dat)
    evs_log = log
    for i in range(len(evs_dat)):
        # log file saves up to 6 decimal points, datfile is full double precision binary
        assert evs_log[i] == pytest.approx(evs_log[i], 1e-6)

def test_hdr_bools(header):
    assert header['eigenfuncs_written']
    assert header['matrices_written']

def test_params(ds):
    expected_params = {'k2': 0, 'k3': 1, 'beta': 0.25,
                       'cte_rho0': 1.0, 'cte_B02': 0.0, 'cte_B03': 1.0}
    assert ds.parameters == expected_params

def test_density(ds):
    rho = ds.equilibria['rho0']
    drho = ds.equilibria['drho0']
    assert (rho == 1).all()
    assert (drho == 0).all()

def test_temperature(ds):
    T0 = ds.equilibria['T0']
    dT0 = ds.equilibria['dT0']
    assert (T0 == 0.125).all()
    assert (dT0 == 0).all()

def test_bfield(ds):
    B02 = ds.equilibria['B02']
    dB02 = ds.equilibria['dB02']
    B03 = ds.equilibria['B03']
    dB03 = ds.equilibria['dB03']
    B0 = ds.equilibria['B0']
    assert (B02 == 0).all()
    assert (dB02 == 0).all()
    assert (B03 == 1).all()
    assert (dB03 == 0).all()
    assert (B0 == 1).all()

def test_velocity(ds):
    v02 = ds.equilibria['v02']
    dv02 = ds.equilibria['dv02']
    v03 = ds.equilibria['v03']
    dv03 = ds.equilibria['dv03']
    assert (v02 == 0).all()
    assert (dv02 == 0).all()
    assert (v03 == 0).all()
    assert (dv03 == 0).all()

def test_cooling_terms(ds):
    dLdT = ds.equilibria['dLdT']
    dLdrho = ds.equilibria['dLdrho']
    assert dLdT.all() == 0
    assert dLdrho.all() == 0

def test_conduction_terms(ds):
    kappa_para = ds.equilibria['kappa_para']
    kappa_perp = ds.equilibria['kappa_perp']
    assert (kappa_para == 0).all()
    assert (kappa_perp == 0).all()

def test_resistivity(ds):
    eta = ds.equilibria['eta']
    detadT = ds.equilibria['detadT']
    assert (eta == 1e-3).all()
    assert (detadT == 0).all()

def test_gravity(ds):
    grav = ds.equilibria['grav']
    assert (grav == 0).all()

def test_speeds(ds):
    cs = ds.get_sound_speed()
    vA = ds.get_alfven_speed()
    cs_expected = np.sqrt(5/3 * 0.125)
    assert (cs == cs_expected).all()
    assert (vA == 1).all()

def test_k0_squared(ds):
    k0_sq = ds.get_k0_squared()
    assert k0_sq == 1

