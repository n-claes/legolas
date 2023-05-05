import numpy as np
import pylbo
import pytest

from .viz_modes import ModeVizTest


class Slice3D(ModeVizTest):
    @property
    def xlabel(self):
        pass

    @property
    def ylabel(self):
        pass

    @property
    def zlabel(self):
        pass

    @property
    def u2vals(self):
        pass

    @property
    def u3vals(self):
        pass

    @property
    def omega(self):
        pass

    @property
    def background(self):
        return False

    @pytest.fixture(scope="function")  # function scope to avoid caching
    def view(self, ds):
        p = pylbo.plot_3d_slice(
            ds,
            omega=self.omega,
            ef_name="rho",
            u2=self.u2vals,
            u3=self.u3vals,
            time=0,
            slicing_axis=self.slicing_axis,
            add_background=self.background,
        )
        p.draw()
        return p

    def test_cbar_lims(self, view, mode_solution):
        assert self.cbar_matches(view, mode_solution)

    def test_labels(self, view):
        assert view.ax.get_xlabel() == self.xlabel
        assert view.ax.get_ylabel() == self.ylabel
        assert view.ax.get_zlabel() == self.zlabel
        assert view.cbar.ax.get_xlabel() == "Re($\\rho$)"

    def test_animation(self, view, tmpdir, mode_solution):
        view.create_animation(
            times=np.arange(5), filename=tmpdir / "test_3d.mp4", fps=1
        )
        assert view.update_colorbar is True
        assert np.allclose(view.solutions, mode_solution)


class TestSliceZ_3DCart(Slice3D):
    filename = "slice_3d_z_cart_rho.npy"
    omega = 1.19029 + 3.75969j
    slicing_axis = "z"
    u2vals = np.linspace(0, 1, 25)
    u3vals = np.linspace(0, 1, 5)
    xlabel = "x"
    ylabel = "y"
    zlabel = "z"

    @pytest.fixture(scope="class")
    def ds(self, ds_v121_rti_khi):
        return ds_v121_rti_khi

    def test_animation_cbar_lock(self, view, tmpdir, mode_solution):
        view.update_colorbar = False
        view.create_animation(
            times=np.arange(5), filename=tmpdir / "test_3d_cbar_lock.mp4", fps=1
        )
        assert view.update_colorbar is False
        assert self.cbar_matches(view, mode_solution)
        assert np.allclose(view.solutions, mode_solution)


class TestSliceZ_3DCyl(Slice3D):
    filename = "slice_3d_z_cyl_rho.npy"
    omega = 0.01746995 + 0.02195201j
    slicing_axis = "z"
    u2vals = np.linspace(0, 2 * np.pi, 25)
    u3vals = np.arange(5)
    xlabel = "x"
    ylabel = "y"
    zlabel = "z"

    @pytest.fixture(scope="class")
    def ds(self, ds_v121_magth):
        return ds_v121_magth
