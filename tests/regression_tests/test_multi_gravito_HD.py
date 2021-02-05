import pytest
import numpy as np
import pylbo
from regression_tests.suite_utils import FIG_DPI, baseline_dir

NB_RUNS = 24

gravito_hd_setup = {
    "name": "gravito_HD",
    "config": {
        "equilibrium_type": "gravito_acoustic",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": np.linspace(0, np.sqrt(250), NB_RUNS),
            "k3": np.linspace(0, np.sqrt(250), NB_RUNS),
            "cte_p0": 1,
            "g": 0.5,
            "alpha": 20.42,
        },
        "basename_datfile": "gravito_acoustic",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
}
parametrisation = dict(
    argnames="setup",
    argvalues=[gravito_hd_setup],
    ids=[gravito_hd_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
@pytest.mark.mpl_image_compare(
    baseline_dir=str(baseline_dir),
    filename="gravito_acoustig.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_gravito_hd(series_test):
    xdata = series_test.get_k0_squared()
    y_scaling = 1 / series_test.get_sound_speed(which_values="average") ** 2
    p = pylbo.plot_spectrum_multi(
        series_test, xdata=xdata, use_squared_omega=True, markersize=3
    )
    p.set_y_scaling(y_scaling)
    p.ax.set_xlim(0, 550)
    p.ax.set_ylim(0, 550)
    return p.fig
