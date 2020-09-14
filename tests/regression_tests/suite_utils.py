import pytest
import numpy as np
import matplotlib.pyplot as plt
import pylbo

# relative tolerance, this is percentage-based (1% here)
REL_TOL = 1e-2
# absolute tolerance
ABS_TOL = 1e-10
# pixel tolerance
PIX_TOL = 1
output = (pylbo.LEGOLAS_DIR / "tests/regression_tests/results").resolve()
answers = (pylbo.LEGOLAS_DIR / "tests/regression_tests/answers").resolve()


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


def compare_eigenvalues(values_test, values_answer, ds_name):
    # first sort values: first on real part, then on imaginary part
    values_test = np.sort_complex(values_test)
    values_answer = np.sort_complex(values_answer)

    fig, ax = plt.subplots(1, figsize=(12, 8))
    fail_idxs = []
    for i, val_test in enumerate(values_test):
        # calculate distance from this point to all points
        distances = np.sqrt(
            (values_answer.real - val_test.real) ** 2
            + (values_answer.imag - val_test.imag) ** 2
        )
        # find point with closest distance
        match = values_answer[distances.argmin()]
        # answer + test in pixels
        x_answ, y_answ = ax.transData.transform((match.real, match.imag))
        x_test, y_test = ax.transData.transform((val_test.real, val_test.imag))
        # get distance between both in pixels
        dist_pix = np.sqrt((x_answ - x_test) ** 2 + (y_answ - y_test) ** 2)
        try:
            # distance in pixels should lie below tolerance
            assert dist_pix < PIX_TOL
        except AssertionError:
            fail_idxs.append(i)
    ax.clear()

    # this should be empty if all values match. If not, raise error and save log
    if fail_idxs:
        for i in fail_idxs:
            val_test = values_test[i]
            val_answ = values_answer[i]
            print(f"FAIL: {val_test} >> {val_answ} (test >> answer)")
            ax.plot(val_test.real, val_test.imag, ".g", markersize=3)
            ax.plot(val_answ.real, val_answ.imag, "xr", markersize=3)
        ax.set_title(ds_name, fontsize=15)
        ax.set_xlabel(r"Re($\omega$)", fontsize=15)
        ax.set_ylabel(r"Im($\omega$)", fontsize=15)
        ax.tick_params(which="both", labelsize=15)
        ax.legend(["test results", "answer results"], loc="best", fontsize=13)
        print(f">>> FAILED EIGENVALUES: {len(fail_idxs)}/{len(values_answer)}")

        filename = (output / f"FAILED_{ds_name}.png").resolve()
        fig.savefig(filename, dpi=400)
        plt.close(fig)
        raise AssertionError
    plt.close(fig)


def compare_eigenfunctions(values_test, values_answer, use_abs=True):
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
    for val_test_arr, val_answ_arr in zip(
        [vals_test_real, vals_test_imag], [vals_answ_real, vals_answ_imag]
    ):
        for val_test, val_answ in zip(val_test_arr, val_answ_arr):
            assert val_test == pytest.approx(val_answ, rel=REL_TOL, abs=ABS_TOL)
