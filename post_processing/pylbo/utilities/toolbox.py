import copy
import numpy as np
from pylbo.utilities.logger import pylboLogger


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
    start : :class:`int`
        The current index in `iterable`, incremented with `step`.
    itr : :class:`~typing.Iterable`
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
