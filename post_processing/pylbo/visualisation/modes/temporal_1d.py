import numpy as np
from pylbo.visualisation.modes.mode_data import ModeVisualisationData
from pylbo.visualisation.modes.mode_figure import ModeFigure2D


class TemporalEvolutionPlot1D(ModeFigure2D):
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
        super().__init__(figsize, data, **kwargs)
        self._u2 = u2
        self._u3 = u3
        self._time = time
        # transpose here so ef_data[:, i] gives eigenfunction at time i
        ef_data = np.broadcast_to(
            self.data.eigenfunction, shape=(len(time), len(self.data.eigenfunction))
        ).transpose()
        time_data = np.broadcast_to(time, shape=ef_data.shape)
        self.set_plot_data(
            u1_data=self.data.ds.ef_grid,
            u2_data=u2,
            u3_data=u3,
            ef_data=ef_data,
            t_data=time_data,
        )
        self._kwargs = kwargs

    def add_mode_solution(self) -> None:
        """Adds the eigenmode solution to the figure."""
        im = self.ax.imshow(
            self.solutions.transpose(),  # transpose due to the nature of imshow
            extent=[
                np.min(self.u1_data),
                np.max(self.u1_data),
                np.min(self.t_data),
                np.max(self.t_data),
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
            The label for the y-axis.
        """
        return "time"
