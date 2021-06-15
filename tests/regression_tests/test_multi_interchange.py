import pytest
import numpy as np
import pylbo
from regression_tests.suite_utils import FIG_DPI

baseline_dir = ""

NB_RUNS = 24

interchange_setup = {
    "name": "interchange",
    "config": {
        "equilibrium_type": "interchange_modes",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": np.pi * np.sin(np.linspace(0, np.pi, NB_RUNS)),
            "k3": np.pi * np.cos(np.linspace(0, np.pi, NB_RUNS)),
            "g": 0.5,
            "cte_p0": 0.25,
            "lambda": 0.3,
            "alpha": 20.0,
        },
        "external_gravity": True,
        "basename_datfile": "interchange",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
}
parametrisation = dict(
    argnames="setup",
    argvalues=[interchange_setup],
    ids=[interchange_setup["name"]],
)


@pytest.mark.skip(reason="needs update")
@pytest.mark.parametrize(**parametrisation)
@pytest.mark.mpl_image_compare(
    baseline_dir=str(baseline_dir),
    filename="interchange_modes.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_interchange_flow(series_test):
    xdata = np.linspace(0, np.pi, len(series_test)) / np.pi
    y_scaling = 1 / series_test.get_alfven_speed(which_values="average") ** 2
    p = pylbo.plot_spectrum_multi(
        series_test, xdata=xdata, use_squared_omega=True, markersize=3
    )
    p.set_y_scaling(y_scaling)
    p.ax.set_xlim(-0.01, 1.01)
    p.ax.set_ylim(-4.1, 14.4)
    p.ax.set_xticks(np.arange(0, 1.2, 0.2))
    p.ax.set_yticks(np.arange(-4, 16, 2))
    return p.fig
