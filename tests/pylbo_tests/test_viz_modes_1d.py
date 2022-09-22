import numpy as np
import pylbo
import pytest

from .viz_modes import ModeVizTest


class TestTemporal1DCart1Mode(ModeVizTest):
    filename = "temporal_1d_cart_1mode_rho.npy"

    @pytest.fixture(scope="class")
    def view(self, ds_v121_rti_khi):
        p = pylbo.plot_1d_temporal_evolution(
            ds_v121_rti_khi,
            omega=1.19029136 + 3.75969744j,
            ef_name="rho",
            u2=1,
            u3=1,
            time=np.linspace(0, 2, 100),
        )
        p.draw()
        return p

    def test_1d_mode_cbar_lims(self, view, mode_solution):
        assert self.cbar_matches(view, mode_solution)

    def test_1d_ylabel(self, view):
        assert view.ax.get_ylabel() == "time"

    def test_1d_cbar_label(self, view):
        assert view.cbar.ax.get_ylabel() == "Re($\\rho$)"


class TestTemporal1dCyl1Mode(TestTemporal1DCart1Mode):
    filename = "temporal_1d_cyl_1mode_rho.npy"

    @pytest.fixture(scope="class")
    def view(self, ds_v121_magth):
        p = pylbo.plot_1d_temporal_evolution(
            ds_v121_magth,
            omega=0.01746995 + 0.02195201j,
            ef_name="rho",
            u2=1,
            u3=1,
            time=np.linspace(0, 100, 100),
        )
        p.draw()
        return p


class TestTemporal1dCart3Modes(TestTemporal1DCart1Mode):
    filename = "temporal_1d_cart_3modes_rho.npy"

    @pytest.fixture(scope="class")
    def view(self, ds_v121_rti_khi):
        p = pylbo.plot_1d_temporal_evolution(
            ds_v121_rti_khi,
            omega=[1.19029 + 3.75969j, 0.65179 + 1.32900j, 0.16193 + 0.45623j],
            ef_name="rho",
            u2=1,
            u3=1,
            time=np.linspace(0, 2, 100),
        )
        p.draw()
        return p
