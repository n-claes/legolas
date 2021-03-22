import numpy as np
from copy import copy
from matplotlib import colors
from pylbo.visualisation.figure_manager import FigureWindow
from pylbo.visualisation.eigenfunctions import EigenfunctionHandler
from pylbo.visualisation.continua import ContinuaHandler
from pylbo.visualisation.legend_interface import LegendHandler
from pylbo.utilities.toolbox import transform_to_numpy, add_pickradius_to_item


class SpectrumFigure(FigureWindow):
    """Superclass of both single and multispectra."""

    def __init__(self, figure_type, figsize, custom_figure):
        super().__init__(
            figure_type=figure_type, figsize=figsize, custom_figure=custom_figure
        )
        self._c_handler = None
        self._ef_handler = None
        self._ef_ax = None

        self.plot_props = None
        self.marker = None
        self.color = None
        self.markersize = None
        self.alpha = None

    def draw(self):
        """
        Draws everything, checks for continua/eigenfunctions. Overridden by subclasses.
        """
        super().draw()
        self._add_spectrum()
        if self._c_handler is not None:
            self.add_continua(self._c_handler.interactive)
        if self._ef_handler is not None:
            self.add_eigenfunctions()

    def _set_plot_properties(self, properties):
        """
        Sets all relevant plot properties.

        Parameters
        ----------
        properties : dict
            Dictionary containing the usual matplotlib properties (marker, color,
            markersize, alpha, etc.)
        """
        plot_props = copy(properties)
        self.marker = plot_props.pop("marker", ".")
        self.color = plot_props.pop("color", "blue")
        self.markersize = plot_props.pop("markersize", 6)
        self.alpha = plot_props.pop("alpha", 0.8)
        self.plot_props = plot_props

    def _add_spectrum(self):
        """Adds the spectrum, is overridden in subclasses."""
        pass

    def add_continua(self, interactive):
        """
        Adds the continua to the plot, overridden by subclasses.

        Parameters
        ----------
        interactive : bool
            If `True`, makes the legend interactive.
        """
        if interactive:
            super()._enable_interactive_legend(self._c_handler)

    def add_eigenfunctions(self):
        """
        Method to add eigenfunctions, should be partially overridden by subclass and
        then called through `super()` to connect the figure events.
        """
        callback_kinds = ("pick_event", "key_press_event")
        callback_methods = (
            self._ef_handler.on_point_pick,
            self._ef_handler.on_key_press,
        )
        for callback_kind, callback_method in zip(callback_kinds, callback_methods):
            callback_id = self.fig.canvas.mpl_connect(callback_kind, callback_method)
            self._mpl_callbacks.append(
                {"cid": callback_id, "kind": callback_kind, "method": callback_method}
            )


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
            figure_type="single-spectrum", figsize=figsize, custom_figure=custom_figure
        )
        self.dataset = dataset
        super()._set_plot_properties(kwargs)

        self.w_real = self.dataset.eigenvalues.real
        self.w_imag = self.dataset.eigenvalues.imag
        self._add_spectrum()

    def _add_spectrum(self):
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
            min_value = np.min(continuum)
            max_value = np.max(continuum)
            # check if continua are complex
            if np.all(np.iscomplex(continuum)):
                item = self.ax.scatter(
                    continuum.real * self.x_scaling,
                    continuum.imag * self.y_scaling,
                    marker=".",
                    color=color,
                    linewidth=self.markersize,
                    alpha=self.alpha / 1.5,
                    label=key,
                )
            elif self._c_handler.check_if_collapsed(continuum=continuum):
                item = self.ax.scatter(
                    min_value * self.x_scaling,
                    0,
                    marker=self._c_handler.marker,
                    s=self._c_handler.markersize,
                    c=color,
                    alpha=self._c_handler.alpha_point,
                    label=key,
                )
            else:
                props = dict(
                    facecolor=colors.to_rgba(color, self._c_handler.alpha_region),
                    edgecolor=colors.to_rgba(color, self._c_handler.alpha_point),
                    label=key,
                )
                if key == "thermal":
                    item = self.ax.axhspan(
                        min_value.imag * self.y_scaling,
                        max_value * self.y_scaling,
                        **props,
                    )
                else:
                    item = self.ax.axvspan(
                        min_value * self.x_scaling, max_value * self.x_scaling, **props
                    )
            self._c_handler.add(item)
        self._c_handler.legend = self.ax.legend(**self._c_handler.legend_properties)
        super().add_continua(self._c_handler.interactive)
        return self._c_handler

    def add_eigenfunctions(self):
        """Adds the eigenfunctions to the plot, sets the eigenfunction handler."""
        if self._ef_ax is None:
            self._ef_ax = super()._add_subplot_axes(self.ax, loc="right")
        if self._ef_handler is None:
            self._ef_handler = EigenfunctionHandler(self.dataset, self._ef_ax)
        # connect everything
        super().add_eigenfunctions()

    def _ensure_smooth_continuum(self, continuum):
        # TODO: this method should split the continuum into multiple parts, such that
        #       regions that lie far apart are not connected by a line.
        pass

    @property
    def ef_ax(self):
        """Property, returns the eigenfunction axes."""
        return self._ef_ax

    @property
    def ef_handler(self):
        """Property, returns the eigenfunction handler."""
        return self._ef_handler


class MultiSpectrumPlot(SpectrumFigure):
    """
    Subclass that draws the multispectra.

    Parameters
    ----------
    dataseries : ~pylbo.data_containers.LegolasDataSeries
        The dataseries that should be used.
    xdata : str, list, numpy.ndarray
        Data to use for the horizontal axis. This can either be a key from the
        parameters dictionary, or a list/numpy array containing actual data.
    use_squared_omega : bool
        If `True`, this will square the eigenvalues when they are plotted on the
        vertical axis. If `False` (default), either the real or imaginary part of the
        eigenvalues will be plotted depending on the value of `use_real_parts`.
    use_real_parts : bool
        If `True` (default), this will plot the real part of the eigenvalues on the
        vertical axis. If `False` the imaginary part will be used.
    figsize : tuple
        Optional figure size like the usual matplotlib (x, x) size.
    custom_figure : tuple
        Optional, in the form (fig, ax). If supplied no new figure will be created
        but this one will be used instead. `fig` refers to the matplotlib figure and
        `ax` to a (single) axes instance, meaning that you can pass a subplot as well.
    """

    def __init__(
        self,
        dataseries,
        xdata,
        use_squared_omega,
        use_real_parts,
        figsize,
        custom_figure,
        **kwargs,
    ):
        super().__init__(
            figure_type="multi-spectrum", figsize=figsize, custom_figure=custom_figure
        )
        self.dataseries = dataseries
        self.use_squared_omega = use_squared_omega
        self._w_pow = 1
        if self.use_squared_omega:
            self._w_pow = 2
        self.use_real_parts = use_real_parts
        self.xdata = self._validate_xdata(xdata)
        self.ydata = self._get_ydata()
        self.x_scaling = np.ones_like(self.dataseries, dtype=float)
        self.y_scaling = np.ones_like(self.dataseries, dtype=float)
        super()._set_plot_properties(kwargs)
        self._add_spectrum()

    def _validate_xdata(self, xdata):
        """
        Validates the xdata passed, does typechecking and necessary casting.
        If a string is passed, this will request the proper values based on the
        parameters.

        Parameters
        ----------
        xdata : str, list, numpy.ndarray
            The xdata used as x values on the spectrum plot.

        Returns
        -------
        xdata_values : numpy.ndarray
            The xdata values of proper length and casted to a Numpy array.
        """
        if isinstance(xdata, str):
            if self.dataseries.parameters.get(xdata, None) is None:
                raise ValueError(
                    f"Provided key xdata='{xdata}' is not in parameters: \n"
                    f"{self.dataseries.parameters.keys()}"
                )
            xdata_values = self.dataseries.parameters[xdata]
        elif isinstance(xdata, (list, np.ndarray)):
            xdata_values = transform_to_numpy(xdata)
            if len(xdata_values) != len(self.dataseries):
                raise ValueError(
                    f"Lengts of xdata do not match: "
                    f"{len(xdata_values)} vs {len(self.dataseries)}"
                )
        else:
            raise TypeError(
                f"xdata should be a string, list or numpy array but got {type(xdata)}."
            )
        return xdata_values

    def _get_ydata(self):
        """
        Gets the y data based on the value of :attr:`use_squared_omega`.

        Returns
        -------
        ydata_values : numpy.ndarray
            The y data values, either the real or imaginary parts based on
            :attr:`use_real_parts`.
        """
        ydata_values = np.array(
            [ds.eigenvalues ** self._w_pow for ds in self.dataseries]
        )
        if self.use_real_parts:
            ydata = ydata_values.real
        else:
            ydata = ydata_values.imag
        # omit zeros from data
        ydata[np.where(np.isclose(ydata, 0, atol=1e-12))] = np.nan
        return ydata

    def set_x_scaling(self, x_scaling):
        """
        Sets the x scaling, properly adjusted to the dataseries length.

        Parameters
        ----------
        x_scaling : int, float, complex, numpy.ndarray
            Values to use for the x-scaling.
        """
        if isinstance(x_scaling, (int, float, complex)):
            x_scaling = np.ones_like(self.dataseries, dtype=float) * x_scaling
        super().set_x_scaling(x_scaling)

    def set_y_scaling(self, y_scaling):
        """
        Sets the y scaling, properly adjusted to the dataseries length.

        Parameters
        ----------
        y_scaling : int, float, complex, numpy.ndarray
            Values to use for the y-scaling.
        """
        if isinstance(y_scaling, (int, float, complex)):
            y_scaling = np.ones_like(self.dataseries, dtype=float) * y_scaling
        super().set_y_scaling(y_scaling)

    def _add_spectrum(self):
        """
        Draw method, creates the spectrum.
        """
        for i, ds in enumerate(self.dataseries):
            (spectrum_point,) = self.ax.plot(
                self.xdata[i]
                * np.ones_like(self.ydata[i], dtype=float)
                * self.x_scaling[i],
                self.ydata[i] * self.y_scaling[i],
                marker=self.marker,
                color=self.color,
                markersize=self.markersize,
                alpha=self.alpha,
                linestyle="None",
                **self.plot_props,
            )
            add_pickradius_to_item(item=spectrum_point, pickradius=10)
            # set dataset associated with this line of points
            setattr(spectrum_point, "dataset", ds)
        self.ax.axhline(y=0, linestyle="dotted", color="grey", alpha=0.3)
        self.ax.axvline(x=0, linestyle="dotted", color="grey", alpha=0.3)

    def add_continua(self, interactive=True):
        """
        Adds the continua to the plot, either interactive or not.

        Parameters
        ----------
        interactive : bool
            If `True`, makes the legend interactive.
        """
        if self._c_handler is None:
            self._c_handler = ContinuaHandler(interactive=interactive)

        for key, color in zip(
            self._c_handler.continua_names, self._c_handler.continua_colors
        ):
            # we skip duplicates if eigenvalues are squared
            if self.use_squared_omega:
                if key in ("slow-", "alfven-"):
                    continue
            # retrieve continuum, calculate region boundaries
            continuum = self.dataseries.continua[key] ** self._w_pow
            if self.use_real_parts:
                continuum = continuum.real
            else:
                continuum = continuum.imag
            min_values = np.array([np.min(c_ds) for c_ds in continuum])
            max_values = np.array([np.max(c_ds) for c_ds in continuum])
            # skip if continua are all zero
            if all(np.isclose(min_values, 0)) and all(np.isclose(max_values, 0)):
                continue
            # when continua are collapsed then min = max and we draw a line instead
            if all(np.isclose(abs(max_values - min_values), 0)):
                (item,) = self.ax.plot(
                    self.xdata * self.x_scaling,
                    min_values * self.y_scaling,
                    color=color,
                    alpha=self._c_handler.alpha_point,
                    label=key,
                )
            else:
                item = self.ax.fill_between(
                    self.xdata * self.x_scaling,
                    min_values * self.y_scaling,
                    max_values * self.y_scaling,
                    facecolor=colors.to_rgba(color, self._c_handler.alpha_region),
                    edgecolor=colors.to_rgba(color, self._c_handler.alpha_point),
                    linewidth=self._c_handler.linewidth,
                    label=key,
                )
            self._c_handler.add(item)
        self._c_handler.legend = self.ax.legend(**self._c_handler.legend_properties)
        super().add_continua(self._c_handler.interactive)
        return self._c_handler

    def add_eigenfunctions(self):
        """Adds the eigenfunctions to the current figure."""
        if self._ef_ax is None:
            self._ef_ax = super()._add_subplot_axes(self.ax, loc="right")
        if self._ef_handler is None:
            self._ef_handler = EigenfunctionHandler(self.dataseries, self._ef_ax)
        # connect everything
        super().add_eigenfunctions()


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
            figure_type="merged-spectrum", figsize=figsize, custom_figure=custom_figure
        )
        self.data = data
        self.leg_handle = LegendHandler(interactive)
        super()._set_plot_properties(kwargs)
        self._use_legend = legend
        self._single_color = False
        if isinstance(kwargs.get("color", None), str):
            self._single_color = True
            # if everything is 1 color no use for a legend
            self._use_legend = False
        self._add_spectrum()

        if self._use_legend and interactive:
            self._enable_interactive_legend(self.leg_handle)

    def _add_spectrum(self):
        """Adds the spectrum to the plot, makes the points pickable."""
        color = None
        if self._single_color:
            color = self.color
        for ds in self.data:
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
        self.fig.tight_layout()

    def add_eigenfunctions(self):
        """Adds the eigenfunctions to the current figure."""
        if self._ef_ax is None:
            self._ef_ax = super()._add_subplot_axes(self.ax, loc="right")
        if self._ef_handler is None:
            self._ef_handler = EigenfunctionHandler(self.data, self._ef_ax)
            # connect everything
        super().add_eigenfunctions()

    def add_continua(self, interactive=True):
        raise NotImplementedError("Continua are not supported for this type of figure.")


class SpectrumComparisonPlot(SpectrumFigure):
    """
    Subclass to compare two spectra.

    Parameters
    ----------
    ds1 : ~pylbo.data_containers.LegolasDataSet
        First dataset, will be placed on the left side.
    ds2 : ~pylbo.data_containers.LegolasDataSet
        Second dataset for comparison, will be placed on the right side.
    figsize : tuple
        Figure size used when creating a window, analogous to matplotlib.
    custom_figure : tuple
        The custom figure to use in the form (fig, axes).
    lock_zoom : bool
        If `True`, locks the zoom for both spectrum plots.

    Attributes
    ----------
    ax2 : ~matplotlib.axes.Axes
        The (top)right axes object.
    panel1 : ~pylbo.visualisation.spectra.SingleSpectrumPlot
        The spectrum instance associated with the left side.
    panel2 : ~pylbo.visualisation.spectra.SingleSpectrumPlot
        The spectrum instance associated with the right side.
    """

    def __init__(self, ds1, ds2, figsize, custom_figure, lock_zoom, **kwargs):
        super().__init__(
            figure_type="compare-spectra",
            figsize=figsize,
            custom_figure=custom_figure,
        )
        super()._set_plot_properties(kwargs)
        share = None
        if lock_zoom:
            share = "all"
        self.ax2 = super()._add_subplot_axes(self.ax, "right", share=share)
        # both panels are essentially single spectra, so create two instances and
        # link that figure with this one
        self.panel1 = SingleSpectrumPlot(
            ds1, figsize=figsize, custom_figure=(self.fig, self.ax)
        )
        self.panel1.ax.set_title(ds1.datfile.stem)
        self.panel2 = SingleSpectrumPlot(
            ds2, figsize=figsize, custom_figure=(self.fig, self.ax2)
        )
        self.panel2.ax.set_title(ds2.datfile.stem)
        self._axes_set = False

    def _use_custom_axes(self):
        """Splits the original 1x2 plot into a 2x2 plot."""
        if self._axes_set:
            return
        self.panel1.ax.change_geometry(2, 2, 1)
        self.panel1._ef_ax = self.panel1.fig.add_subplot(2, 2, 3)
        self.panel2.ax.change_geometry(2, 2, 2)
        self.panel2._ef_ax = self.panel2.fig.add_subplot(2, 2, 4)
        self._axes_set = True

    def add_eigenfunctions(self):
        """Adds the eigenfunctions for both datasets and merges the mpl callbacks."""
        self._use_custom_axes()
        for panel in [self.panel1, self.panel2]:
            panel.add_eigenfunctions()
            panel.disconnect_callbacks()
            # merge callbacks
            self._mpl_callbacks.extend(panel._mpl_callbacks)
            # add dedicated attribute to prevent mpl from double triggering events
            setattr(panel.ef_handler, "associated_ds_ax", panel.ax)

    def add_continua(self, interactive=True):
        """Adds the continua for both datasets and merges the mpl callbacks."""
        for panel in (self.panel1, self.panel2):
            panel.add_continua(interactive)
            panel.disconnect_callbacks()
            self._mpl_callbacks.extend(panel._mpl_callbacks)
