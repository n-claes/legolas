import functools
import time

import matplotlib.lines as mpl_lines
import numpy as np
from pylbo._version import _mpl_version
from pylbo.utilities.logger import pylboLogger


def timethis(func):
    @functools.wraps(func)
    def _time_method(*args, **kwargs):
        t0 = time.perf_counter()
        try:
            return func(*args, **kwargs)
        finally:
            pylboLogger.debug(
                f"{func.__name__} took {time.perf_counter() - t0} seconds to execute"
            )

    return _time_method


def get_axis_geometry(ax):
    """
    Retrieves the geometry of a given matplotlib axis.

    Parameters
    ----------
    ax : ~matplotlib.axes.Axes
        The axis to retrieve the geometry from.

    Returns
    -------
    tuple
        The geometry of the given matplotlib axis.
    """
    if _mpl_version >= "3.4":
        axis_geometry = ax.get_subplotspec().get_geometry()[0:3]
    else:
        # this is 1-based indexing by default, use 0-based here for consistency
        # with subplotspec in matplotlib 3.4+
        axis_geometry = transform_to_numpy(ax.get_geometry())
        axis_geometry[-1] -= 1
        axis_geometry = tuple(axis_geometry)
    return axis_geometry


def add_pickradius_to_item(item, pickradius):
    """
    Makes a matplotlib artist pickable and adds a pickradius.
    We have to handle this separately, because for line2D items the method
    :meth:`~matplotlib.axes.Axes.set_picker` is deprecated from version 3.3 onwards.

    Parameters
    ----------
    item : ~matplotlib.artist.Artist
        The artist which will be made pickable
    pickradius : int, float
        Sets the pickradius, which determines if something is "on" the picked point.
    """
    # set_picker is deprecated for line2D from matplotlib 3.3 onwards
    if isinstance(item, mpl_lines.Line2D) and _mpl_version >= "3.3":
        item.set_picker(True)
        item.pickradius = pickradius
    else:
        item.set_picker(pickradius)


def custom_enumerate(iterable, start=0, step=1):
    """
    Does a custom enumeration with a given stepsize.

    Parameters
    ----------
    iterable : ~typing.Iterable
        The iterable to iterate over.
    start : int
        The starting value for enumerate.
    step : int
        The stepsize between enumerate values.

    Yields
    ------
    start : int
        The current index in `iterable`, incremented with `step`.
    itr : ~typing.Iterable
        The corresponding entry of `iterable`.
    """
    for itr in iterable:
        yield start, itr
        start += step


def transform_to_list(obj: any) -> list:
    """
    Transforms a given input argument `obj` to a list. If `obj`
    is a Numpy-array or tuple, a cast to `list()` is invoked.

    Parameters
    ----------
    obj : any
        The object to transform.

    Returns
    -------
    list
        The object converted to a list.
    """
    if obj is None:
        return [obj]
    elif isinstance(obj, (tuple, np.ndarray)):
        return list(obj)
    elif isinstance(obj, list):
        return obj
    return [obj]


def transform_to_numpy(obj: any) -> np.ndarray:
    """
    Transforms a given input argument `obj` to a numpy array.

    Parameters
    ----------
    obj : any
        The object to transform.

    Returns
    -------
    numpy.ndarray
        The object transformed to a numpy array.

    """
    if obj is None:
        return np.asarray([obj])
    elif isinstance(obj, (tuple, list)):
        return np.asarray(obj)
    elif isinstance(obj, np.ndarray):
        return obj
    return np.asarray([obj])


def solve_cubic_exact(a, b, c, d):
    """
    Solves a given cubic polynomial of the form
    :math:`ax^3 + bx^2 + cx + d = 0` using the analytical cubic root formula
    instead of the general `numpy.roots` routine.
    From `StackOverflow <https://math.stackexchange.com/questions
    15865why-not-write-the-solutions-of-a-cubic-this-way/18873#18873/>`_.

    Parameters
    ----------
    a : int, float, complex
        Cubic coefficient.
    b : int, float, complex
        Quadratic coefficient.
    c : int, float, complex
        Linear coefficient.
    d : int, float, complex
        Constant term

    Returns
    -------
    roots : np.ndarray(ndim=3, dtype=complex)
        The three roots of the cubic polynomial as a Numpy array.
    """

    if a == 0:
        raise ValueError("cubic coefficient may not be zero")
    p = b / a
    q = c / a
    r = d / a
    Aterm = (
        -2 * p**3
        + 9 * p * q
        - 27 * r
        + 3
        * np.sqrt(3)
        * np.sqrt(
            -(p**2) * q**2
            + 4 * q**3
            + 4 * p**3 * r
            - 18 * p * q * r
            + 27 * r**2
        )
    ) ** (1 / 3) / (3 * 2 ** (1 / 3))
    Bterm = (-(p**2) + 3 * q) / (9 * Aterm)
    cterm_min = (-1 - np.sqrt(3) * 1j) / 2
    cterm_pos = (-1 + np.sqrt(3) * 1j) / 2
    x1 = -p / 3 + Aterm - Bterm
    x2 = -p / 3 + cterm_min * Aterm - cterm_pos * Bterm
    x3 = -p / 3 + cterm_pos * Aterm - cterm_min * Bterm
    return np.array([x1, x2, x3], dtype=complex)
