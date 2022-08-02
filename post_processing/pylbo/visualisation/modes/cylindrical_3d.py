import numpy as np
from pylbo.visualisation.modes.cartesian_3d import CartesianSlicePlot3D
from pylbo.visualisation.modes.mode_data import ModeVisualisationData


class CylindricalSlicePlot3D(CartesianSlicePlot3D):
    """
    Class for handling cylindrical 3D plots of the eigenmode solution.

    Parameters
    ----------
    data : ModeVisualisationData
        The data for the visualisation.
    u2 : np.ndarray
        The :math:`\\theta` coordinate of the eigenmode solution.
    u3 : np.ndarray
        The :math:`z` coordinate of the eigenmode solution.
    time : float
        The time at which the eigenmode solution is calculated.
    slicing_axis : str
        The axis along which the eigenmode solution is sliced.
    figsize : tuple[int, int]
        The size of the figure.
    vmin : float
        The minimum value of the colourbar. If None, the minimum value of the
        solution is used.
    vmax : float
        The maximum value of the colourbar. If None, the maximum value of the
        solution is used.
    """

    def __init__(
        self,
        data: ModeVisualisationData,
        u2: np.ndarray,
        u3: np.ndarray,
        time: float,
        slicing_axis: str,
        figsize: tuple[int, int],
        vmin: float = None,
        vmax: float = None,
        **kwargs,
    ) -> None:
        if figsize is None:
            figsize = (8, 8)
        if slicing_axis == "theta":
            raise NotImplementedError(
                "3D slicing is not implemented for a theta slice."
            )
        super().__init__(data, u2, u3, time, slicing_axis, figsize, **kwargs)

        self.vmin = np.min(self._solutions) if vmin is None else vmin
        self.vmax = np.max(self._solutions) if vmax is None else vmax

    def set_plot_arrays(self) -> None:
        ef = self.data.eigenfunction
        rgrid, thetagrid = np.meshgrid(self.data.ds.ef_grid, self._u2, indexing="ij")
        self.ef_data = np.broadcast_to(ef, shape=(len(self._u2), len(ef))).transpose()
        self.u1_data = rgrid * np.cos(thetagrid)
        self.u2_data = rgrid * np.sin(thetagrid)
        self.u3_data = self._u3
        self.time_data = self._time

    def draw_solution(self) -> None:
        super().draw_solution()
        xmax = np.max(self._u1)
        self.ax.set_xlim(-xmax, xmax)
        self.ax.set_ylim(-xmax, xmax)

    def get_view_xlabel(self) -> str:
        return "x"

    def get_view_ylabel(self) -> str:
        return "y"
