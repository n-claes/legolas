from pylbo._version import __version__
from pylbo.utilities import logger
from pylbo.utilities.logger import set_loglevel, disable_logging
from pylbo.file_handler import load, load_series, load_logfile

logger.init_logger()
