import warnings
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.figure import Figure
from pylbo.visualisation.modes.mode_data import ModeVisualisationData
from pylbo.visualisation.modes.mode_figure import ModeFigure2D
from pylbo.visualisation.utils import add_axis_label


class SpatialCartesianPlot2D(ModeFigure2D):
    """
    Class for handling Cartesian 2D plots of the eigenmode solutions.

    Parameters
    ----------
    data : ModeVisualisationData
        The data for the visualisation.
    u2 : float or ndarray
        The :math:`y`  coordinate of the eigenmode solution.
    u3 : float or ndarray
        The :math:`z`  coordinate of the eigenmode solution.
    time : float
        The time at which the eigenmode solution is calculated.
    slicing_axis : str
        The axis along which the eigenmode solution is sliced.
    figsize : tuple[int, int]
        The size of the figure.
    **kwargs
        Additional keyword arguments to be passed to :meth:`matplotlib.pyplot.imshow`.
    """

    def __init__(
        self,
        data: ModeVisualisationData,
        u2: Union[float, np.ndarray],
        u3: Union[float, np.ndarray],
        time: float,
        slicing_axis: str,
        figsize: tuple[int, int],
        **kwargs,
    ) -> None:
        for i in "123":
            _axis = getattr(data.ds, f"u{i}_str")
            setattr(self, f"_u{i}axis", _axis.replace("$", "").replace("\\", ""))
        self.slicing_axis = self._validate_slicing_axis(
            slicing_axis, allowed_axes=[self._u2axis, self._u3axis]
        )
        self._u2 = self._validate_u2(u2, slicing_axis, axis=self._u2axis)
        self._u3 = self._validate_u3(u3, slicing_axis, axis=self._u3axis)
        self._time = time
        super().__init__(figsize, data)

        self._set_plot_data_slice()
        self._kwargs = kwargs

    def _set_plot_data_slice(self) -> None:
        """Sets the plot data for a slice along a given axis."""
        ef = self.data.eigenfunction
        axis = self.slicing_axis
        coord = self._u2 if axis == self._u3axis else self._u3
        ef_2d = np.broadcast_to(ef, shape=(len(coord), len(ef))).transpose()
        x_2d, coord_2d = np.meshgrid(self.data.ds.ef_grid, coord, indexing="ij")
        self.set_plot_data(
            u1_data=x_2d,
            u2_data=coord_2d if axis == self._u3axis else self._u2,
            u3_data=coord_2d if axis == self._u2axis else self._u3,
            ef_data=ef_2d,
            t_data=self._time,
        )

    def add_u2u3_txt(self, ax, **kwargs) -> None:
        """
        Creates a textbox on the figure with the value of the :math:`u_2-u_3`
        coordinates.

        Parameters
        ----------
        ax : ~matplotlib.axes.Axes
            The axes to use for the textbox.
        **kwargs
            Additional keyword arguments to pass to :meth:`add_axis_label`.
        """
        if self.slicing_axis == self._u3axis:
            txt = rf"{self.data.ds.u3_str} = {self._u3}"
        else:
            txt = rf"{self.data.ds.u2_str} = {self._u2}"
        txt = rf"{txt} | t = {self._time}"
        self.u2u3_txt = add_axis_label(ax, txt, **kwargs)

    def add_mode_solution(self) -> None:
        """Adds the eigenmode solution to the figure"""
        if self.slicing_axis == self._u3axis:
            extent_vertical = (np.min(self._u2), np.max(self._u2))
        else:
            extent_vertical = (np.min(self._u3), np.max(self._u3))
        im = self.ax.imshow(
            self.solutions.transpose(),
            extent=[
                np.min(self.data.ds.ef_grid),
                np.max(self.data.ds.ef_grid),
                *extent_vertical,
            ],
            aspect="auto",
            origin="lower",
            **self._kwargs,
        )
        self.cbar = self.fig.colorbar(im, cax=self.cbar_ax)

    def get_view_ylabel(self) -> str:
        """
        Returns
        -------
        str
            The label for the y-axis on the bottom panel of the figure.
        """
        return (
            self.data.ds.u2_str
            if self.slicing_axis == self._u3axis
            else self.data.ds.u3_str
        )


class SpatialCylindricalPlot2D(SpatialCartesianPlot2D):
    """
    Class for handling cylindrical 2D plots of the eigenmode solutions.

    Parameters
    ----------
    data : ModeVisualisationData
        The data for the visualisation.
    u2 : float or ndarray
        The :math:`\\theta`  coordinate of the eigenmode solution.
    u3 : float or ndarray
        The :math:`z`  coordinate of the eigenmode solution.
    time : float
        The time at which the eigenmode solution is calculated.
    slicing_axis : str
        The axis along which the eigenmode solution is sliced.
    figsize : tuple[int, int]
        The size of the figure.
    **kwargs
        Additional keyword arguments to be passed to
        :meth:`matplotlib.pyplot.pcolormesh`.
    """

    def __init__(
        self,
        data: ModeVisualisationData,
        u2: Union[float, np.ndarray],
        u3: Union[float, np.ndarray],
        time: float,
        slicing_axis: str,
        figsize: tuple[int, int],
        polar: bool,
        **kwargs,
    ) -> None:
        self._use_polar_axes = polar
        super().__init__(data, u2, u3, time, slicing_axis, figsize, **kwargs)

    def _set_plot_data_slice(self) -> None:
        """Sets the plot data for a slice along a given axis."""
        if self.slicing_axis == self._u2axis:
            return super()._set_plot_data_slice()
        ef = self.data.eigenfunction
        thetas = self._u2
        ef_2d = np.broadcast_to(ef, shape=(len(thetas), len(ef))).transpose()
        r_2d, theta_2d = np.meshgrid(self.data.ds.ef_grid, thetas, indexing="ij")
        self.set_plot_data(
            u1_data=r_2d,
            u2_data=theta_2d,
            u3_data=self._u3,
            ef_data=ef_2d,
            t_data=self._time,
        )

    def add_mode_solution(self) -> None:
        """Adds the eigenmode solution to the figure"""
        if self.slicing_axis == self._u2axis:
            return super().add_mode_solution()
        x_2d = self.u1_data * np.cos(self.u2_data)
        y_2d = self.u1_data * np.sin(self.u2_data)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            if not self._use_polar_axes:
                im = self.ax.pcolormesh(x_2d, y_2d, self.solutions, **self._kwargs)
            else:
                im = self.ax.pcolormesh(
                    self.u2_data, self.u1_data, self.solutions, **self._kwargs
                )
        self.cbar = self.fig.colorbar(
            ScalarMappable(norm=im.norm, cmap=im.cmap), cax=self.cbar_ax
        )

    def add_eigenfunction(self) -> None:
        """Adds the eigenfunction to the figure."""
        super().add_eigenfunction()
        self.axes["eigfunc"].set_xlabel(self.data.ds.u1_str)

    def get_view_xlabel(self) -> str:
        """
        Returns
        -------
        str
            The label for the x-axis on the bottom panel of the figure.
        """
        if self.slicing_axis == self._u3axis:
            return "x"
        return super().get_view_xlabel()

    def get_view_ylabel(self) -> str:
        """
        Returns
        -------
        str
            The label for the y-axis on the bottom panel of the figure.
        """
        if self._use_polar_axes:
            return ""
        if self.slicing_axis == self._u3axis:
            return "y"
        return super().get_view_ylabel()

    def _create_figure_layout(self, figsize: tuple[int, int]) -> tuple[Figure, dict]:
        """
        Overloads the superclass method for figure layout creation.

        Parameters
        ----------
        figsize : tuple[int, int]
            The size of the figure.

        Returns
        -------
        fig : ~matplotlib.figure.Figure
            The figure to use for the visualisation.
        axes : dict
            The axes to use for the visualisation.
        """
        if self.slicing_axis == self._u2axis:
            return super()._create_figure_layout(figsize)
        fig = plt.figure(figsize=figsize)
        polar = self._use_polar_axes
        ax1 = fig.add_axes([0.1, 0.7, 0.8, 0.2])
        ax2 = fig.add_axes([0.25, 0.1, 0.5, 0.5], aspect="equal", polar=polar)
        if polar:
            self._cbar_hspace = 0.05
        return fig, {"eigfunc": ax1, "view": ax2}
