from multiprocessing.sharedctypes import Value

import numpy as np
import pylbo
import pytest

from .viz_modes import ModeVizTest


class SliceSetup(ModeVizTest):
    @property
    def xlabel(self):
        pass

    @property
    def ylabel(self):
        pass

    @property
    def u2vals(self):
        pass

    @property
    def u3vals(self):
        pass

    @pytest.fixture(scope="class")
    def view(self, ds_v121_rti_khi):
        p = pylbo.plot_2d_slice(
            ds_v121_rti_khi,
            omega=1.19029 + 3.75969j,
            ef_name="rho",
            u2=self.u2vals,
            u3=self.u3vals,
            time=0,
            slicing_axis=self.slicing_axis,
        )
        p.draw()
        return p

    def test_cbar_lims(self, view, mode_solution):
        actual = (view.cbar.vmin, view.cbar.vmax)
        expected = (mode_solution.min(), mode_solution.max())
        assert np.allclose(actual, expected)

    def test_labels(self, view):
        assert view.ax.get_xlabel() == self.xlabel
        assert view.ax.get_ylabel() == self.ylabel
        assert view.cbar.ax.get_ylabel() == "Re($\\rho$)"


class Slice2D(SliceSetup):
    def test_invalid_slicing_axis(self, ds_v121_rti_khi):
        with pytest.raises(ValueError):
            pylbo.plot_2d_slice(
                ds_v121_rti_khi,
                omega=1.19029 + 3.75969j,
                ef_name="rho",
                u2=np.linspace(0, 1, 10),
                u3=1,
                time=0,
                slicing_axis="x",
            )


class TestSliceZ_2DCart(Slice2D):
    filename = "slice_2d_z_cart_rho.npy"
    slicing_axis = "z"
    u2vals = np.linspace(0, 2, 50)
    u3vals = 0
    xlabel = "x"
    ylabel = "y"

    def test_invalid_coords(self, ds_v121_rti_khi):
        with pytest.raises(ValueError):
            pylbo.plot_2d_slice(
                ds_v121_rti_khi,
                omega=1.19029 + 3.75969j,
                ef_name="rho",
                u2=0,
                u3=np.linspace(0, 1, 10),
                time=0,
                slicing_axis="z",
            )


class TestSliceZ_2DCyl(Slice2D):
    filename = "slice_2d_z_cyl_rho.npy"
    slicing_axis = "z"
    u2vals = np.linspace(0, 2 * np.pi, 50)
    u3vals = 1
    xlabel = "x"
    ylabel = "y"


class TestSliceY_2DCart(Slice2D):
    filename = "slice_2d_y_cart_rho.npy"
    slicing_axis = "y"
    u2vals = 1
    u3vals = np.linspace(0, 2, 50)
    xlabel = "x"
    ylabel = "z"

    def test_invalid_coords(self, ds_v121_rti_khi):
        with pytest.raises(ValueError):
            pylbo.plot_2d_slice(
                ds_v121_rti_khi,
                omega=1.19029 + 3.75969j,
                ef_name="rho",
                u2=np.linspace(0, 1, 10),
                u3=0,
                time=0,
                slicing_axis="y",
            )
