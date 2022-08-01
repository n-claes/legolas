import warnings
from typing import Union

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.figure import Figure
from pylbo.visualisation.modes.cartesian_2d import CartesianSlicePlot2D
from pylbo.visualisation.modes.mode_data import ModeVisualisationData


class CylindricalSlicePlot2D(CartesianSlicePlot2D):
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
    polar : bool
        Whether to use polar coordinates for the plot.
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

    def set_plot_arrays(self) -> None:
        if self.slicing_axis == self._u2axis:
            return super().set_plot_arrays()
        ef = self.data.eigenfunction
        thetas = self._u2
        self.ef_data = np.broadcast_to(ef, shape=(len(thetas), len(ef))).transpose()
        r_2d, theta_2d = np.meshgrid(self.data.ds.ef_grid, thetas, indexing="ij")

        self.u1_data = r_2d
        self.u2_data = theta_2d
        self.u3_data = self._u3
        self.time_data = self._time

    def draw_solution(self) -> None:
        """Adds the eigenmode solution to the figure"""
        if self.slicing_axis == self._u2axis:
            return super().draw_solution()
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

    def draw_eigenfunction(self) -> None:
        """Adds the eigenfunction to the figure."""
        super().draw_eigenfunction()
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
