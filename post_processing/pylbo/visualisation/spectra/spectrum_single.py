import numpy as np
from pylbo.utilities.toolbox import add_pickradius_to_item
from pylbo.visualisation.continua import ContinuaHandler
from pylbo.visualisation.eigenfunctions.derived_eigfunc_handler import (
    DerivedEigenfunctionHandler,
)
from pylbo.visualisation.eigenfunctions.eigfunc_handler import EigenfunctionHandler
from pylbo.visualisation.spectra.spectrum_figure import SpectrumFigure


class SingleSpectrumPlot(SpectrumFigure):
    """
    Creates a plot of a single spectrum based on a given dataset.

    Parameters
    ----------
    dataset : ~pylbo.data_containers.LegolasDataSet
        The dataset used to create the spectrum.
    figsize : tuple
        Figure size used when creating a window, analogous to matplotlib.
    custom_figure : tuple
        The custom figure to use in the form (fig, axes).

    Attributes
    ----------
    dataset : ~pylbo.data_containers.LegolasDataSet
        The dataset passed as parameter
    w_real : numpy.ndarray(dtype=float, ndim=1)
        Real part of the eigenvalues as a numpy array.
    w_imag : numpy.ndarray(dtype=float, ndim=1)
        Imaginary part of the eigenvalues as a numpy array.
    marker : ~matplotlib.markers
        The marker used to draw the points.
    markersize : int, float
        Size of the marker.
    alpha : int, float
        Alpha value of the points.
    """

    def __init__(self, dataset, figsize, custom_figure, **kwargs):
        super().__init__(
            custom_figure=custom_figure, figlabel="single-spectrum", figsize=figsize
        )
        self.dataset = dataset
        super()._set_plot_properties(kwargs)

        self.w_real = self.dataset.eigenvalues.real
        self.w_imag = self.dataset.eigenvalues.imag

    def add_spectrum(self):
        """Adds the spectrum to the plot, makes the points pickable."""
        (spectrum_point,) = self.ax.plot(
            self.w_real * self.x_scaling,
            self.w_imag * self.y_scaling,
            marker=self.marker,
            color=self.color,
            markersize=self.markersize,
            alpha=self.alpha,
            linestyle="None",
            **self.plot_props,
        )
        # set dataset associated with this line of points
        setattr(spectrum_point, "dataset", self.dataset)
        add_pickradius_to_item(item=spectrum_point, pickradius=10)
        self.ax.axhline(y=0, linestyle="dotted", color="grey", alpha=0.3)
        self.ax.axvline(x=0, linestyle="dotted", color="grey", alpha=0.3)
        self.ax.set_xlabel(r"Re($\omega$)")
        self.ax.set_ylabel(r"Im($\omega$)")
        self.ax.set_title(self.dataset.eq_type)

    def add_continua(self, interactive=True):
        """
        Adds the continua to the spectrum.

        Parameters
        ----------
        interactive : bool
            If `True`, makes the legend pickable.

        Returns
        -------
        c_handler : ~pylbo.continua.ContinuaHandler
            The legendhandler used to plot the continua.
        """
        if self._c_handler is None:
            self._c_handler = ContinuaHandler(interactive=interactive)

        for key, color in zip(
            self._c_handler.continua_names, self._c_handler.continua_colors
        ):
            continuum = self.dataset.continua[key]
            if self._c_handler.check_if_all_zero(continuum=continuum):
                continue
            realvals = continuum.real * self.x_scaling
            imagvals = continuum.imag * self.y_scaling
            if self._c_handler.check_if_collapsed(continuum=continuum):
                realvals = np.min(realvals)
                imagvals = 0
            item = self.ax.scatter(
                realvals,
                imagvals,
                marker=self._c_handler.marker,
                linewidth=self._c_handler.markersize,
                c=color,
                alpha=self._c_handler.alpha_point,
                label=key,
            )
            self._c_handler.add(item)
        self._c_handler.legend = self.ax.legend(**self._c_handler.legend_properties)
        if interactive:
            super().make_legend_interactive(self._c_handler)

    def add_eigenfunctions(self):
        """Adds the eigenfunctions to the plot, sets the eigenfunction handler."""
        if self._ef_ax is None:
            self._ef_ax = super().add_subplot_axes(self.ax, loc="right")
        if self._ef_handler is None:
            self._ef_handler = EigenfunctionHandler(self.dataset, self._ef_ax, self.ax)
        super().add_eigenfunction_interface(efhandler=self._ef_handler)

    def add_derived_eigenfunctions(self):
        """
        Adds the derived eigenfunctions to the plot, sets the eigenfunction handler.
        """
        if self._def_ax is None:
            self._def_ax = super().add_subplot_axes(self.ax, loc="right")
        if self._def_handler is None:
            self._def_handler = DerivedEigenfunctionHandler(
                self.dataset, self._def_ax, self.ax
            )
        super().add_eigenfunction_interface(efhandler=self._def_handler)
