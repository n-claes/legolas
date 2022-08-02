from typing import Union

import numpy as np
from pylbo.visualisation.modes.mode_data import ModeVisualisationData
from pylbo.visualisation.modes.mode_figure import ModeFigure
from pylbo.visualisation.utils import add_axis_label


class CartesianSlicePlot2D(ModeFigure):
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
        self._u1 = data.ds.ef_grid
        self._u2 = self._validate_u2(u2, slicing_axis, axis=self._u2axis)
        self._u3 = self._validate_u3(u3, slicing_axis, axis=self._u3axis)
        self._time = self._check_if_number(time, "time")
        self._kwargs = kwargs
        super().__init__(figsize, data)

    def _validate_slicing_axis(self, slicing_axis: str, allowed_axes: list[str]) -> str:
        """
        Validates the slicing axis.

        Parameters
        ----------
        slicing_axis : str
            The slicing axis.
        allowed_axes : list[str]
            The list of allowed axes.

        Returns
        -------
        str
            The validated slicing axis.
        """
        if slicing_axis not in allowed_axes:
            raise ValueError(f"Slicing axis must be one of {allowed_axes}.")
        return slicing_axis

    def _validate_u2(self, u2: float, slicing_axis: str, axis: str) -> float:
        """
        Validates the combination of u2 and slicing axis.

        Parameters
        ----------
        u2 : float
            The :math:`u_2` coordinate.
        slicing_axis : str
            The slicing axis.
        axis : str
            The coordinate axis corresponding to :math:`u_2`.

        Returns
        -------
        float
            The validated :math:`u_2` coordinate.
        """
        if slicing_axis == axis and not isinstance(u2, (int, np.integer, float)):
            raise ValueError(f"u2 must be a number for slicing axis '{axis}'.")
        return u2

    def _validate_u3(self, u3: float, slicing_axis: str, axis: str) -> float:
        """
        Validates the combination of u3 and slicing axis.

        Parameters
        ----------
        u3 : float
            The :math:`u_3` coordinate.
        slicining_axis : str
            The slicing axis.
        axis : str
            The coordinate axis corresponding to :math:`u_3`.

        Returns
        -------
        float
            The validated :math:`u_3` coordinate.
        """
        if slicing_axis == axis and not isinstance(u3, (int, np.integer, float)):
            raise ValueError(f"u3 must be a number for slicing axis '{axis}'.")
        return u3

    def set_plot_arrays(self) -> None:
        ef = self.data.eigenfunction
        axis = self.slicing_axis
        coord = self._u2 if axis == self._u3axis else self._u3
        self.ef_data = np.broadcast_to(ef, shape=(len(coord), len(ef))).transpose()
        x_2d, coord_2d = np.meshgrid(self.data.ds.ef_grid, coord, indexing="ij")

        self.u1_data = x_2d
        self.u2_data = coord_2d if axis == self._u3axis else self._u2
        self.u3_data = coord_2d if axis == self._u2axis else self._u3
        self.time_data = self._time

    def add_u2u3_txt(self, ax, **kwargs) -> None:
        if self.slicing_axis == self._u3axis:
            txt = rf"{self.data.ds.u3_str} = {self._u3}"
        else:
            txt = rf"{self.data.ds.u2_str} = {self._u2}"
        txt = rf"{txt} | t = {self._time}"
        self.u2u3_txt = add_axis_label(ax, txt, **kwargs)

    def draw_solution(self) -> None:
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
        return (
            self.data.ds.u2_str
            if self.slicing_axis == self._u3axis
            else self.data.ds.u3_str
        )
