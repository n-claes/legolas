from typing import Union

import numpy as np
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
        super().__init__(figsize, data)

        self.slicing_axis = self._validate_slicing_axis(
            slicing_axis, allowed_axes=["z", "y"]
        )
        self._u2 = self._validate_u2(u2, slicing_axis, coord_axis="y")
        self._u3 = self._validate_u3(u3, slicing_axis, coord_axis="z")
        self._time = time
        self._set_plot_data_slice()
        self._kwargs = kwargs

    def _set_plot_data_slice(self) -> None:
        """
        Sets the plot data for a slice along 'z' (xy-plane) or along 'y' (xz-plane).
        """
        ef = self.data.eigenfunction
        axis = self.slicing_axis
        coord = self._u2 if axis == "z" else self._u3
        ef_2d = np.broadcast_to(ef, shape=(len(coord), len(ef))).transpose()
        _, coord_2d = np.meshgrid(self.data.ds.ef_grid, coord, indexing="ij")

        self.set_plot_data(
            u1_data=self.data.ds.ef_grid,
            u2_data=coord_2d if axis == "z" else self._u2,
            u3_data=coord_2d if axis == "y" else self._u3,
            ef_data=ef_2d,
            t_data=self._time,
        )

    def add_u2u3_txt(self, ax, **kwargs) -> None:
        if self.slicing_axis == "z":
            txt = rf"{self.data.ds.u3_str} = {self._u3}"
        else:
            txt = rf"{self.data.ds.u2_str} = {self._u2}"
        txt = rf"{txt} | t = {self._time}"
        self.u2u3_txt = add_axis_label(ax, txt, **kwargs)

    def add_mode_solution(self) -> None:
        """Adds the eigenmode solution to the figure"""
        if self.slicing_axis == "z":
            extent_vertical = (np.min(self._u2), np.max(self._u2))
        else:
            extent_vertical = (np.min(self._u3), np.max(self._u3))
        im = self.ax.imshow(
            self.solutions.transpose(),
            extent=[
                np.min(self.u1_data),
                np.max(self.u1_data),
                *extent_vertical,
            ],
            aspect="auto",
            origin="lower",
            **self._kwargs,
        )
        self.cbar = self.fig.colorbar(im, cax=self.cbar_ax)
        self.cbar.set_label(rf"{self.data._ef_name_latex}")
        self.ax.set_xlabel(self.data.ds.u1_str)
        self.ax.set_ylabel(
            self.data.ds.u2_str if self.slicing_axis == "z" else self.data.ds.u3_str
        )
