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
        super().__init__(figsize)
        self.data = data

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

    def draw(self) -> None:
        """
        Draws everything on the figure: eigenfunction, eigenmode solution,
        and the various text boxes.
        """
        self.add_eigenfunction()
        self.add_mode_solution()

        self.add_omega_txt(self.axes["eigfunc"], loc="top left", outside=True)
        self.add_u2u3_txt(self.axes["eigfunc"], loc="top right", outside=True)
        self.add_k2k3_txt(self.ax, loc="bottom left", color="white", alpha=0.5)

    def add_eigenfunction(self) -> None:
        """Adds the eigenfunction to the figure."""
        ax = self.axes["eigfunc"]
        ef = getattr(self.data.eigenfunction, self.data.part_name)
        ax.plot(self.u1_data, ef, lw=2)
        ax.axvline(x=0, color="grey", ls="--", lw=1)
        ax.set_xlim(np.min(self.u1_data), np.max(self.u1_data))
        ax.set_title(rf"Temporal evolution of {self.data._ef_name_latex}")
        ax.set_ylabel(self.data._ef_name_latex)

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
        self.cbar.set_label(rf"{self.data._ef_name_latex}")
        self.ax.set_xlabel(self.data.ds.u1_str)
        self.ax.set_ylabel("time")
