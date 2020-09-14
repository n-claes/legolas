import pytest
import matplotlib.pyplot as plt
from pathlib import Path
import pylbo

pylbo.set_loglevel("warning")


@pytest.fixture(scope="module", autouse=True)
def ds():
    file = Path("utility_files/v1_datfile_efs.dat").resolve()
    return pylbo.load(file)


def test_passedlist(ds):
    with pytest.raises(TypeError):
        pylbo.SingleSpectrum([ds, ds])


def test_plot_spectrum(ds):
    s = pylbo.SingleSpectrum(ds)
    s.plot_spectrum()
    plt.close(s.fig)


def test_plot_spectrum_passfig(ds):
    fig, ax = plt.subplots(1)
    s = pylbo.SingleSpectrum(ds, fig=fig, ax=ax)
    s.plot_spectrum()
    plt.close(s.fig)


def test_plot_continua(ds):
    s = pylbo.SingleSpectrum(ds)
    s.plot_continua()
    plt.close(s.fig)


def test_plot_equilibria(ds):
    s = pylbo.SingleSpectrum(ds)
    s.plot_equilibria()
    plt.close(s.fig)


def test_plot_efs_nomerge(ds):
    s = pylbo.SingleSpectrum(ds)
    s.plot_eigenfunctions(merge_figs=False)
    plt.close(s.fig)


def test_plot_efs_merge(ds):
    s = pylbo.SingleSpectrum(ds)
    s.plot_eigenfunctions(merge_figs=True)
    plt.close(s.fig)


def test_plot_efs_nocont(ds):
    s = pylbo.SingleSpectrum(ds)
    s.plot_eigenfunctions(merge_figs=True, annotate_continua=False)
    plt.close(s.fig)


def test_inconsistent_run(ds):
    from pylbo.utilities.exceptions import InconsistentMultirunFile

    file2 = Path("utility_files/v1_resistive.dat").resolve()
    ds2 = pylbo.load(file2)
    ms = pylbo.MultiSpectrum([ds, ds2])
    with pytest.raises(InconsistentMultirunFile):
        ms.plot_precoded_run()

        
def test_noprecoded_run(ds):
    ms = pylbo.MultiSpectrum([ds, ds])
    with pytest.raises(ValueError):
        ms.plot_precoded_run()
