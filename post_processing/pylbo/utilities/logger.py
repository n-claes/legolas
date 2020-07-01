import logging
import sys

pylboLogger = logging.getLogger('pylbo')

def init_logger():
    sh = logging.StreamHandler(stream=sys.stdout)
    fmt = logging.Formatter('%(asctime)s - |%(levelname)-8s| %(message)s')
    sh.setFormatter(fmt)
    pylboLogger.addHandler(hdlr=sh)
    set_loglevel("info")

def set_loglevel(level):
    if isinstance(level, str):
        level = level.upper()
    pylboLogger.setLevel(level)

def disable_logging():
    for hdlr in pylboLogger.handlers:
        pylboLogger.removeHandler(hdlr)
    pylboLogger.addHandler(logging.NullHandler())


