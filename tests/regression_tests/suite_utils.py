import pytest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# relative tolerance, this is percentage-based (1% here)
REL_TOL = 1e-2
# absolute tolerance
ABS_TOL = 1e-10
# pixel tolerance
PIX_TOL = 1
# dpi for image-based comparison
FIG_DPI = 200
# eigenfunction names
EF_NAMES = ("rho", "v1", "v2", "v3", "T", "a1", "a2", "a3")

SAVEFIG_KWARGS = {"dpi": FIG_DPI, "transparent": True}
RMS_TOLERANCE = 2

regression_test_dir = Path(__file__).parent
output = (regression_test_dir / "results").resolve()
answers = (regression_test_dir / "answers").resolve()

if not output.is_dir():
    output.mkdir()
if not answers.is_dir():
    answers.mkdir()


def get_filepaths(basename):
    filename_dat = "{}.dat".format(basename)
    filename_log = "{}.log".format(basename)
    datfile = (output / filename_dat).resolve()
    logfile = (output / filename_log).resolve()
    return datfile, logfile


def get_answer_filepaths(basename):
    datfile, logfile = get_filepaths(basename)
    answer_datfile = (answers / "answer_{}".format(datfile.name)).resolve()
    answer_logfile = (answers / "answer_{}".format(logfile.name)).resolve()
    return answer_datfile, answer_logfile


def get_image_filename(name, limits_dict):
    xlims = limits_dict["xlims"]
    ylims = limits_dict["ylims"]
    figname = f"{name}_Re{xlims[0]}-{xlims[1]}_Im{ylims[0]}-{ylims[1]}"
    return figname
