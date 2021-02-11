import pytest
import numpy as np
import pylbo
from regression_tests.suite_utils import FIG_DPI, baseline_dir

NB_RUNS = 24

constant_current_setup = {
    "name": "constant_current",
    "config": {
        "equilibrium_type": "constant_current_tokamak",
        "gridpoints": 51,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": -2.0,
            "k3": 0.2,
            "j0": (2.0 * 0.2) / np.linspace(1.9, 2.1, NB_RUNS),
            "cte_rho0": 1.0,
            "cte_B03": 1.0,
        },
        "basename_datfile": "constant_current",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
}
parametrisation = dict(
    argnames="setup",
    argvalues=[constant_current_setup],
    ids=[constant_current_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
@pytest.mark.mpl_image_compare(
    baseline_dir=str(baseline_dir),
    filename="constant_current.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_constant_current(series_test):
    xdata = 2 * series_test.parameters.get("k3") / series_test.parameters.get("j0")
    p = pylbo.plot_spectrum_multi(
        series_test, xdata=xdata, use_squared_omega=True, markersize=3
    )
    p.ax.set_yscale("symlog", linthresh=1e-8)
    p.ax.set_xlim(1.88, 2.12)
    p.ax.set_ylim(-1e-3, 1e6)
    return p.fig
