from pylbo.visualisation.figure_manager import FigureWindow
from pylbo.utilities.logger import pylboLogger
from pylbo.visualisation.continua import ContinuaRegions


class SingleSpectrumPlot(FigureWindow):
    def __init__(self, dataset, figsize, custom_figure, **kwargs):
        super().__init__(
            figure_type="single-spectrum",
            figsize=figsize,
            custom_figure=custom_figure
        )
        self.dataset = dataset
        self.kwargs = kwargs
        self.xdata = self.dataset.eigenvalues.real
        self.ydata = self.dataset.eigenvalues.imag

    def draw(self):
        pylboLogger.debug("creating single spectrum plot...")
        self.xdata *= self.x_scaling
        self.ydata *= self.y_scaling
        self.ax.plot(
            self.xdata,
            self.ydata,
            marker=self.kwargs.pop("marker", "."),
            color=self.kwargs.pop("color", "blue"),
            markersize=self.kwargs.pop("markersize", 6),
            alpha=self.kwargs.pop("alpha", 0.8),
            linestyle="None",
            **self.kwargs,
        )
        self.ax.axhline(y=0, linestyle="dotted", color="grey", alpha=0.3)
        self.ax.axvline(x=0, linestyle="dotted", color="grey", alpha=0.3)
        self.ax.set_xlabel(r"Re($\omega$)")
        self.ax.set_ylabel(r"Im($\omega$)")
        self.ax.set_title(self.dataset.eq_type)

    def add_continua(self, interactive=True, **kwargs):
        regions = ContinuaRegions(self.dataset, interactive, self, **kwargs)
        pass

    def add_eigenfunctions(self):
        pass
