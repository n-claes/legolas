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

SAVEFIG_KWARGS = {"dpi": FIG_DPI}
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


def compare_eigenfunctions(values_test, values_answer, use_abs=True, relax=False):
    if use_abs:
        vals_test_real, vals_test_imag = (
            np.abs(values_test.real),
            np.abs(values_test.imag),
        )
        vals_answ_real, vals_answ_imag = (
            np.abs(values_answer.real),
            np.abs(values_answer.imag),
        )
    else:
        vals_test_real, vals_test_imag = values_test.real, values_test.imag
        vals_answ_real, vals_answ_imag = values_answer.real, values_answer.imag
    # real parts
    try:
        assert np.all(
            vals_answ_real == pytest.approx(vals_test_real, rel=REL_TOL, abs=ABS_TOL)
        )
    except AssertionError:
        if relax:
            answ_norm = normalise_eigenfunction(vals_answ_real)
            test_norm = normalise_eigenfunction(vals_test_real)
            assert np.all(abs(answ_norm - test_norm) <= 1e-3)
        else:
            raise
    # imaginary parts
    try:
        assert np.all(
            vals_answ_imag == pytest.approx(vals_test_imag, rel=REL_TOL, abs=ABS_TOL)
        )
    except AssertionError:
        if relax:
            answ_norm = normalise_eigenfunction(vals_answ_imag)
            test_norm = normalise_eigenfunction(vals_test_imag)
            assert np.all(abs(answ_norm - test_norm) <= 1e-3)
        else:
            raise


def normalise_eigenfunction(array):
    if np.all(np.isclose(array, 0, atol=ABS_TOL)):
        return array
    # this normalises between 0 and 1
    arr_norm = (array - np.min(array)) / (np.max(array) - np.min(array))
    return arr_norm
