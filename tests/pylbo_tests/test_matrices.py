import pytest
import numpy as np
from pathlib import Path
import pylbo

pylbo.set_loglevel("warning")
FIG_DPI = 300
BASELINE_DIR = str(Path("utility_files/baseline").resolve())
RESULTS_DIR = str(Path("results").resolve())


@pytest.fixture(scope="module", autouse=True)
def ds():
    file = Path("utility_files/v1_datfile_matrices.dat").resolve()
    return pylbo.load(file)


@pytest.fixture(scope="module", autouse=True)
def ds_adiab():
    file = Path("utility_files/v1_adiabatic_car.dat").resolve()
    return pylbo.load(file)


def test_matrixB_notpresent(ds_adiab):
    from pylbo.utilities.exceptions import MatricesNotPresent

    with pytest.raises(MatricesNotPresent):
        ds_adiab.get_matrix_B()


def test_matrixB():
    file = Path("utility_files/v1_datfile_matrices.dat")
    ds = pylbo.load(file)
    rows, cols, vals = ds.get_matrix_B()
    for arr in (rows, cols):
        assert isinstance(arr, np.ndarray)
        assert np.all([isinstance(val, np.integer) for val in arr])
    assert isinstance(vals, np.ndarray)
    assert np.all([isinstance(val, np.float) for val in vals])


def test_matrixA_notpresent(ds_adiab):
    from pylbo.utilities.exceptions import MatricesNotPresent

    with pytest.raises(MatricesNotPresent):
        ds_adiab.get_matrix_A()


def test_matrixA():
    file = Path("utility_files/v1_datfile_matrices.dat")
    ds = pylbo.load(file)
    rows, cols, vals = ds.get_matrix_A()
    for arr in (rows, cols):
        assert isinstance(arr, np.ndarray)
        assert np.all([isinstance(val, np.integer) for val in arr])
    assert isinstance(vals, np.ndarray)
    assert np.all([isinstance(val, np.complex) for val in vals])


@pytest.mark.mpl_image_compare(
    baseline_dir=BASELINE_DIR,
    filename="matrices_adiab.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_matrix_plot(ds):
    fig, axes = pylbo.plot_matrices(ds)
    return fig
