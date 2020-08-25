import copy
import numpy as np
from .defaults import precoded_runs
from .exceptions import UnknownPrecodedRun
from .logger import pylboLogger


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


def transform_to_list(obj):
    """
    Transforms a given input argument `obj` to a list. If `obj`
    is a Numpy-array, :func:`~numpy.ndarray.tolist` is invoked.

    Parameters
    ----------
    obj : object
        The object to transform.

    Returns
    -------
    obj: list
        The object converted to a list.
    """
    if obj is None:
        return [obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, list):
        return obj
    return [obj]


def transform_to_numpy(obj):
    """
    Transforms a given input argument `obj` to a numpy array.

    Parameters
    ----------
    obj : object
        The object to transform.

    Returns
    -------
    obj : numpy.ndarray
        The object transformed to a numpy array.

    """
    if obj is None:
        return np.asarray([obj])
    elif isinstance(obj, np.ndarray):
        return np.array(obj)
    elif isinstance(obj, (tuple, list)):
        return np.asarray(obj)
    return np.asarray([obj])


def get_precoded_run(name):
    """
    Retrieves the configuration dictionary for a precoded run based on
    its name. These are the keys of `precoded_runs`.

    Parameters
    ----------
    name : str
        The name of the precoded run, as a key of `precoded_runs`.

    Returns
    -------
    selected_run : dict
        A :func:`~copy.deepcopy` of the configuration dictionary for the selected
        precoded run, to avoid conflicts with multiple runs in the same script.

    Raises
    ------
    UnknownPrecodedRun
        If `name` is not a key of `precoded_runs`, and is hence unknown.
    """
    try:
        selected_run = precoded_runs[name]
    except KeyError:
        raise UnknownPrecodedRun(name, precoded_runs.keys())
    return copy.deepcopy(selected_run)


def select_precoded_run():
    """
    Prints all available precoded runs such that a selection can be made.

    Returns
    -------
    parfile_dict : dict
        Dictionary containing the chosen configuration.

    """
    av_runs = iter(precoded_runs.keys())
    pylboLogger.info('select one of the following precoded runs:')
    for idx, pr in custom_enumerate(av_runs, start=1, step=3):
        print('{:03d} {:<30}{:03d} {:<30}{:03d} {}'.format(idx, pr, idx + 1, next(av_runs, "<empty>"),
                                                           idx + 2, next(av_runs, '<empty>')))
    pr_in = int(input('\nChoose precoded run: '))
    chosen_pr = list(precoded_runs.keys())[pr_in - 1]
    pylboLogger.info(f'selected run: {chosen_pr}')
    parfile_dict = get_precoded_run(chosen_pr)
    return parfile_dict
