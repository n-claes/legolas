from pylbo.utilities.toolbox import add_pickradius_to_item
from pylbo.visualisation.eigenfunctions.derived_eigfunc_handler import (
    DerivedEigenfunctionHandler,
)
from pylbo.visualisation.eigenfunctions.eigfunc_handler import EigenfunctionHandler
from pylbo.visualisation.legend_handler import LegendHandler
from pylbo.visualisation.spectra.spectrum_figure import SpectrumFigure


class MergedSpectrumPlot(SpectrumFigure):
    """
    Merges the datasets from a given series into a single plot.

    Parameters
    ----------
    data : ~pylbo.data_containers.LegolasDataSeries
        The dataseries which will be merged.
    figsize : tuple
        Figure size used when creating a window, analogous to matplotlib.
    custom_figure : tuple
        The custom figure to use in the form (fig, axes).
    interactive : bool
        If `True` an interactive legend is enabled.
    legend : bool
        If `False` no legend will be drawn.

    Attributes
    ----------
    data : ~pylbo.data_containers.LegolasDataSeries
        The dataseries passed as parameter.
    leg_handle : ~pylbo.visualisation.legend_interface.LegendHandler
        The handler for the legend.
    """

    def __init__(self, data, figsize, custom_figure, interactive, legend, **kwargs):
        super().__init__(
            custom_figure=custom_figure, figlabel="merged_spectrum", figsize=figsize
        )
        self.data = data
        self.leg_handle = LegendHandler(interactive)
        kwargs_plt = dict(kwargs)
        kwargs_plt.pop("colors", None)
        kwargs_plt.pop("color_parameter", None)
        super()._set_plot_properties(kwargs_plt)
        self._use_legend = legend
        self._single_color = False
        if isinstance(kwargs.get("color", None), str):
            self._single_color = True
            # if everything is 1 color no use for a legend
            self._use_legend = False
        if isinstance(kwargs.get("colors", None), dict):
            if isinstance(kwargs.get("color_parameter", None), str):
                self._color_from_parameter = True
                self._colors = kwargs.get("colors", None) 
                self._color_parameter = kwargs.get("color_parameter", None)
        self._interactive = interactive

    def add_spectrum(self):
        """Adds the spectrum to the plot, makes the points pickable."""
        color = None
        if self._single_color:
            color = self.color
        for ds in self.data:
            if self._color_from_parameter:
                color = self._colors[ds.parameters[self._color_parameter]]
            spectrum_point = self.ax.scatter(
                ds.eigenvalues.real * self.x_scaling,
                ds.eigenvalues.imag * self.y_scaling,
                marker=self.marker,
                s=6 * self.markersize,
                c=color,
                alpha=self.alpha,
                label=ds.datfile.stem,
                **self.plot_props,
            )
            setattr(spectrum_point, "dataset", ds)
            add_pickradius_to_item(item=spectrum_point, pickradius=10)
            self.leg_handle.add(spectrum_point)
        self.ax.axhline(y=0, linestyle="dotted", color="grey", alpha=0.3)
        self.ax.axvline(x=0, linestyle="dotted", color="grey", alpha=0.3)
        self.ax.set_xlabel(r"Re($\omega$)")
        self.ax.set_ylabel(r"Im($\omega$)")

        if self._use_legend:
            self.leg_handle.legend = self.ax.legend(loc="best")
            self.leg_handle._make_visible_by_default = True
            if self._interactive:
                super().make_legend_interactive(self.leg_handle)
        self.fig.tight_layout()

    def add_eigenfunctions(self):
        """Adds the eigenfunctions to the current figure."""
        if self._ef_ax is None:
            self._ef_ax = super().add_subplot_axes(self.ax, loc="right")
        if self._ef_handler is None:
            self._ef_handler = EigenfunctionHandler(self.data, self._ef_ax, self.ax)
        super().add_eigenfunction_interface(efhandler=self._ef_handler)

    def add_derived_eigenfunctions(self):
        """Adds the derived eigenfunctions to the current figure."""
        if self._def_ax is None:
            self._def_ax = super().add_subplot_axes(self.ax, loc="right")
        if self._def_handler is None:
            self._def_handler = DerivedEigenfunctionHandler(
                self.data, self._def_ax, self.ax
            )
        super().add_eigenfunction_interface(efhandler=self._def_handler)

    def add_continua(self, interactive=True):
        raise NotImplementedError("Continua are not supported for this type of figure.")
