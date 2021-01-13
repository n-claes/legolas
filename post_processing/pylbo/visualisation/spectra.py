import numpy as np
from pylbo.visualisation.figure_manager import FigureWindow
from pylbo.utilities.logger import pylboLogger
from pylbo.continua import ContinuaHandler


class SingleSpectrumPlot(FigureWindow):
    """
    Creates a plot of a single spectrum based on a given dataset.

    Parameters
    ----------
    dataset : `~pylbo.data_containers.LegolasDataSet` instance
        The dataset used to create the spectrum.
    figsize : tuple
        Figure size used when creating a window, analogous to matplotlib.
    custom_figure : tuple
        The custom figure to use in the form (fig, axes)

    Attributes
    ----------
    dataset : `~pylbo.data_containers.LegolasDataSet`
        The dataset passed as parameter
    xdata : np.ndarray(dtype=float, ndim=1)
        Real part of the eigenvalues as a numpy array.
    ydata : np.ndarray(dtype=float, ndim=1)
        Imaginary part of the eigenvalues as a numpy array.
    """
    def __init__(self, dataset, figsize, custom_figure, **kwargs):
        super().__init__(
            figure_type="single-spectrum", figsize=figsize, custom_figure=custom_figure
        )
        self.dataset = dataset
        self.kwargs = kwargs
        self.xdata = self.dataset.eigenvalues.real
        self.ydata = self.dataset.eigenvalues.imag
        self.draw()

    def draw(self):
        """
        Draw method, creates the spectrum.
        """
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

    def add_continua(self, interactive=True, colors=None, **kwargs):
        """
        Adds the continua to the spectrum.

        Parameters
        ----------
        interactive : bool
            If `True`, makes the legend pickable.
        colors : np.ndarray(dtype=str, ndim=1)
            Array of colors which will be used to plot the continua.

        Returns
        -------
        c_handler : `~pylbo.continua.ContinuaHandler` instance
            The `ContinuaHandler` instance used to plot the continua.
        """
        c_handler = ContinuaHandler(self, interactive)
        c_handler.continua_colors = colors
        c_handler.alpha_point = kwargs.pop("alpha_point", 0.8)
        c_handler.alpha_region = kwargs.pop("alpha_region", 0.2)
        c_handler.alpha_hidden = kwargs.pop("alpha_hidden", 0.05)

        for key, color in zip(c_handler.continua_names, c_handler.continua_colors):
            continuum = self.dataset.continua[key]
            if c_handler.check_if_all_zero(continuum=continuum):
                continue
            min_value = np.min(continuum)
            max_value = np.max(continuum)
            if c_handler.check_if_collapsed(continuum=continuum):
                item = self.ax.scatter(
                    min_value,
                    0,
                    marker=kwargs.pop("marker", "p"),
                    s=kwargs.pop("markersize", 64),
                    c=color,
                    alpha=c_handler.alpha_point,
                    label=key,
                )
            else:
                props = dict(facecolor=color, alpha=c_handler.alpha_region, label=key)
                if key == "thermal":
                    item = self.ax.axhspan(min_value, max_value, **props)
                else:
                    item = self.ax.axvspan(min_value, max_value, **props)
            c_handler.add(item)
        c_handler.legend = self.ax.legend(**kwargs)

        super().add_continua(c_handler, interactive, kwargs.pop("pickradius", 10))
        return c_handler

    def add_eigenfunctions(self):
        pass
