import numpy as np
from pylbo.visualisation.modes.mode_data import ModeVisualisationData
from pylbo.visualisation.modes.mode_figure import ModeFigure


class TemporalEvolutionPlot1D(ModeFigure):
    """
    Main class for temporal evolutions of the eigenfunction.

    Parameters
    ----------
    data : ModeVisualisationData
        Data object containing all data associated with the selected eigenmode.
    u2 : float
        The data for the :math:`u_2` coordinate.
    u3 : float
        The data for the :math:`u_3` coordinate.
    time : np.ndarray
        The data for the time.
    figsize : tuple[int, int]
        The size of the figure.
    **kwargs
        Additional keyword arguments to be passed to :meth:`matplotlib.pyplot.imshow`.
    """

    def __init__(
        self,
        data: ModeVisualisationData,
        u2: float,
        u3: float,
        time: np.ndarray,
        figsize: tuple[int, int],
        **kwargs,
    ) -> None:
        self._u1 = data.ds.ef_grid
        self._u2 = self._check_if_number(u2, "u2")
        self._u3 = self._check_if_number(u3, "u3")
        self._time = self._check_if_array(time, "time")
        self._kwargs = kwargs
        super().__init__(figsize, data)

    def set_plot_arrays(self) -> None:
        # transpose here so ef_data[:, i] gives eigenfunction at time i
        self.ef_data = np.broadcast_to(
            self.data.eigenfunction,
            shape=(len(self._time), len(self.data.eigenfunction)),
        ).transpose()
        self.time_data = np.broadcast_to(self._time, shape=self.ef_data.shape)
        self.u1_data = self.data.ds.ef_grid
        self.u2_data = self._u2
        self.u3_data = self._u3

    def draw_solution(self) -> None:
        self._view = self.ax.imshow(
            self.solutions.transpose(),  # transpose due to the nature of imshow
            extent=[
                np.min(self.u1_data),
                np.max(self.u1_data),
                np.min(self.time_data),
                np.max(self.time_data),
            ],
            aspect="auto",
            origin="lower",
            **self._kwargs,
        )
        self.cbar = self.fig.colorbar(self._view, cax=self.cbar_ax)

    def get_view_ylabel(self) -> str:
        return "time"
