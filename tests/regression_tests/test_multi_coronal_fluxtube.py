import pytest
import numpy as np
import pylbo
from regression_tests.suite_utils import FIG_DPI, baseline_dir

NB_RUNS = 24

coronal_tube_setup = {
    "name": "coronal_flux_tubes",
    "config": {
        "equilibrium_type": "coronal_flux_tube",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": 0.0,
            "k3": np.linspace(0.1, 6.2, NB_RUNS),
            "cte_rho0": 1.0,
            "cte_p0": 1.0,
            "r0": 1.0,
        },
        "basename_datfile": "coronal_tube",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
}
parametrisation = dict(
    argnames="setup",
    argvalues=[coronal_tube_setup],
    ids=[coronal_tube_setup["name"]],
)


@pytest.mark.parametrize(**parametrisation)
@pytest.mark.mpl_image_compare(
    baseline_dir=str(baseline_dir),
    filename="coronal_fluxtubes.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_coronal_fluxtubes(series_test):
    xdata = (
        series_test.parameters.get("k3")
        * coronal_tube_setup["config"]["parameters"]["r0"]
    )
    y_scaling = 1 / (xdata * series_test.get_sound_speed("maximum"))
    p = pylbo.plot_spectrum_multi(
        series_test, xdata=xdata, use_squared_omega=False, markersize=3
    )
    p.set_y_scaling(y_scaling)
    p.ax.set_ylim(0.84, 5.1)
    return p.fig
