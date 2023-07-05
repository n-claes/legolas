import matplotlib.pyplot as plt
import numpy as np
import pytest
from pylbo.utilities import toolbox


def test_geometry_single_figure():
    _, ax = plt.subplots(1)
    assert toolbox.get_axis_geometry(ax) == (1, 1, 0)


def test_geometry_multiple_columns():
    _, axes = plt.subplots(1, 3)
    assert toolbox.get_axis_geometry(axes[1]) == (1, 3, 1)


def test_geometry_multiple_rows():
    _, axes = plt.subplots(4, 1)
    assert toolbox.get_axis_geometry(axes[2]) == (4, 1, 2)


def test_geometry_multiple_rows_and_columns():
    _, axes = plt.subplots(3, 3)
    assert toolbox.get_axis_geometry(axes[0, 1]) == (3, 3, 1)
    assert toolbox.get_axis_geometry(axes[1, 1]) == (3, 3, 4)
    assert toolbox.get_axis_geometry(axes[2, 0]) == (3, 3, 6)


def test_enumerate():
    expected = [0, 3, 6, 9, 12]
    result = [i for i, _ in toolbox.custom_enumerate(np.arange(5), step=3)]
    assert expected == result


def test_float_tolist():
    result = toolbox.transform_to_list(5)
    assert isinstance(result, list)
    assert result == [5]


def test_cmplx_tolist():
    result = toolbox.transform_to_list(1 + 2j)
    assert isinstance(result, list)
    assert result == [1 + 2j]


def test_numpy_tolist():
    result = toolbox.transform_to_list(np.linspace(0, 1, 6))
    assert isinstance(result, list)
    assert result == pytest.approx([0, 0.2, 0.4, 0.6, 0.8, 1])


def test_list_tolist():
    result = toolbox.transform_to_list([1, 2, 3, 4])
    assert isinstance(result, list)
    assert result == [1, 2, 3, 4]


def test_none_tolist():
    result = toolbox.transform_to_list(None)
    assert isinstance(result, list)
    assert result[0] is None


def test_float_tonumpy():
    result = toolbox.transform_to_numpy(8)
    assert isinstance(result, np.ndarray)
    assert result[0] == 8


def test_cmplx_tonumpy():
    result = toolbox.transform_to_numpy(3 - 5j)
    assert isinstance(result, np.ndarray)
    assert result[0] == 3 - 5j


def test_numpy_tonumpy():
    expected = np.linspace(0, 1, 10)
    result = toolbox.transform_to_numpy(expected)
    assert isinstance(result, np.ndarray)
    assert result == pytest.approx(expected)


def test_list_tonumpy():
    result = toolbox.transform_to_numpy([0, 0.25, 0.5, 0.75, 1])
    assert isinstance(result, np.ndarray)
    assert result == pytest.approx(np.linspace(0, 1, 5))


def test_none_tonumpy():
    result = toolbox.transform_to_numpy(None)
    assert isinstance(result, np.ndarray)
    assert result[0] is None


def test_cubic_solver_a_zero():
    with pytest.raises(ValueError):
        toolbox.solve_cubic_exact(a=0, b=1, c=2, d=3)


def test_cubic_solver():
    sols = np.sort_complex(toolbox.solve_cubic_exact(a=2.5, b=-2, c=1, d=7.5))
    assert all(
        np.isclose(
            sols,
            np.array(
                [-1.14371037, 0.97185518 - 1.2955845j, 0.97185518 + 1.2955845j],
                dtype=complex,
            ),
        )
    )
    assert all(np.isclose(sols, np.sort_complex(np.roots([2.5, -2, 1, 7.5]))))


def test_nzeroes_real(ds_v114_subset):
    evs = np.array([6.7, 11, 16, 21, 26, 31], dtype=complex)
    all_efs = ds_v114_subset.get_eigenfunctions(ev_guesses=evs)
    rho_efs = np.array([efs["rho"] for efs in all_efs], dtype=complex)
    expected_zeroes = np.arange(1, 7)
    result = toolbox.count_zeroes(rho_efs, real=True)
    assert len(result) == len(expected_zeroes)
    assert np.all(result == expected_zeroes)


def test_nzeroes_imag(ds_v114_subset):
    evs = np.array([6.7, 11, 16, 21, 26, 31], dtype=complex)
    all_efs = ds_v114_subset.get_eigenfunctions(ev_guesses=evs)
    vx_efs = np.array([efs["v1"] for efs in all_efs], dtype=complex)
    expected_zeroes = np.arange(0, 6)
    result = toolbox.count_zeroes(vx_efs, real=False)
    assert len(result) == len(expected_zeroes)
    assert np.all(result == expected_zeroes)


def test_nzeroes_allzero():
    efs = np.zeros((3, 10), dtype=complex)
    result = toolbox.count_zeroes(efs, real=True)
    expected_zeros = np.zeros(3, dtype=int)
    assert np.all(result == expected_zeros)


def test_resonance_location_exact():
    grid = np.linspace(0, 1, 100)
    continuum = np.linspace(4, 8, 100)
    sigma = continuum[47]
    result = toolbox.find_resonance_location(continuum, grid, sigma)
    assert np.isclose(grid[47], result)


def test_resonance_location_interp():
    grid = np.linspace(0, 1, 11)
    continuum = np.linspace(-5, 5, 11)
    sigma = -1.5
    result = toolbox.find_resonance_location(continuum, grid, sigma)
    assert np.isclose(0.35, result)


def test_resonance_location_none():
    grid = np.linspace(0, 1, 11)
    continuum = np.linspace(-3, 3, 11)
    result = toolbox.find_resonance_location(continuum, grid, sigma=5)
    assert result is None
    result = toolbox.find_resonance_location(continuum, grid, sigma=-4.5)
    assert result is None


def test_resonance_location_continuum_not_monotone():
    continuum = np.array([0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0])
    grid = np.linspace(0, 1, len(continuum))
    sigma = 3.5
    result = toolbox.find_resonance_location(continuum, grid, sigma)
    assert np.isclose(0.35, result)
