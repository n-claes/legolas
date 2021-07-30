import pytest
import numpy as np
from pylbo.utilities import toolbox


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
