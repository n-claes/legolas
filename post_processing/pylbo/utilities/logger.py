import logging
import sys

pylboLogger = logging.getLogger('pylbo')


def init_logger():
    """
    Initialises the pylbo logger.
    """
    sh = logging.StreamHandler(stream=sys.stdout)
    fmt = logging.Formatter('%(asctime)s - |%(levelname)-8s| %(message)s')
    sh.setFormatter(fmt)
    pylboLogger.addHandler(hdlr=sh)
    set_loglevel("info")


def set_loglevel(level):
    """
    Sets the logging level.

    Parameters
    ----------
    level : int, str
        The level for logging. See :class:`~logging.Logger`, or the
        `allowed levels <https://docs.python.org/3/library/logging.html#levels>`_.
        Both the string and integer values can be set, case-insensitive.
    """
    if isinstance(level, str):
        level = level.upper()
    pylboLogger.setLevel(level)


def disable_logging():
    """
    Completely disables logging.
    """
    for hdlr in pylboLogger.handlers:
        pylboLogger.removeHandler(hdlr)
    pylboLogger.addHandler(logging.NullHandler())
