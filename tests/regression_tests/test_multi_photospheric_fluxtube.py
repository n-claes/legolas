import pytest
import numpy as np
import pylbo
from regression_tests.suite_utils import FIG_DPI

baseline_dir = ""

NB_RUNS = 24

photospheric_tube_setup = {
    "name": "photospheric_flux_tubes",
    "config": {
        "equilibrium_type": "photospheric_flux_tube",
        "gridpoints": 31,
        "number_of_runs": NB_RUNS,
        "parameters": {
            "k2": 0.0,
            "k3": np.linspace(0.1, 6.2, NB_RUNS),
            "cte_rho0": 1.0,
            "cte_p0": 1.0,
            "r0": 1.0,
        },
        "basename_datfile": "photospheric_tube",
        "write_eigenfunctions": False,
        "show_results": False,
        "logging_level": 0,
    },
}
parametrisation = dict(
    argnames="setup",
    argvalues=[photospheric_tube_setup],
    ids=[photospheric_tube_setup["name"]],
)


@pytest.mark.skip(reason="needs update")
@pytest.mark.parametrize(**parametrisation)
@pytest.mark.mpl_image_compare(
    baseline_dir=str(baseline_dir),
    filename="photospheric_fluxtubes.png",
    savefig_kwargs={"dpi": FIG_DPI},
)
def test_photospheric_fluxtubes(series_test):
    xdata = (
        series_test.parameters.get("k3")
        * photospheric_tube_setup["config"]["parameters"]["r0"]
    )
    y_scaling = 1 / (xdata * series_test.get_sound_speed("minimum"))
    p = pylbo.plot_spectrum_multi(
        series_test, xdata=xdata, use_squared_omega=False, markersize=3
    )
    p.set_y_scaling(y_scaling)
    p.ax.set_ylim(0.84, 1.55)
    return p.fig
