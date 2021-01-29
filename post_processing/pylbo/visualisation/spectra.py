import numpy as np
from copy import copy
from matplotlib import colors
from pylbo.visualisation.figure_manager import FigureWindow
from pylbo.visualisation.eigenfunctions import EigenfunctionHandler
from pylbo.utilities.toolbox import transform_to_numpy, add_pickradius_to_item
from pylbo.continua import ContinuaHandler


class SpectrumFigure(FigureWindow):
    def __init__(self, figure_type, figsize, custom_figure):
        super().__init__(
            figure_type=figure_type, figsize=figsize, custom_figure=custom_figure
        )
        self._c_handler = None
        self._ef_handler = None
        self._ef_ax = None

    def draw(self):
        super().draw()
        self._add_spectrum()
        if self._c_handler is not None:
            self.add_continua(self._c_handler.interactive)
        if self._ef_handler is not None:
            self.add_eigenfunctions()

    def _add_spectrum(self):
        pass

    def add_continua(self, interactive):
        """
        Method to add continua, should be partially overridden by subclass and then
        called through `super()`. The attribute `self._c_handler` should have
        been set in the subclass.
        This solely makes an existing legend pickable if `interactive=True`.
        """
        if interactive:
            self._c_handler.make_legend_pickable()
            callback_kind = "pick_event"
            callback_method = self._c_handler.on_legend_pick
            callback_id = self.fig.canvas.mpl_connect(callback_kind, callback_method)
            self._mpl_callbacks.append(
                {
                    "cid": callback_id,
                    "kind": callback_kind,
                    "method": callback_method,
                }
            )

    def _add_eigenfunction_axes(self):
        """
        Adds a new axes to the existing plot.
        This sets the attribute `self._ef_handler` and will modify the geometry of
        the existing axes, adding a new subplot in which the eigenfunctions will
        be drawn.
        Eventually `fig.tight_layout()` is called to update the figure's gridspec.
        """
        if self._ef_ax is not None:
            return

        if self.ax.get_geometry() != (1, 1, 1):
            raise ValueError(
                f"Something went wrong when adding the eigenfunctions. Expected "
                f"axes with geometry (1, 1, 1) but got {self.ax.get_geometry()}."
            )
        self.ax.change_geometry(1, 2, 1)
        self._ef_ax = self.fig.add_subplot(122)
        self.fig.set_size_inches(16, 8)
        # this will update the figure's gridspec
        self.fig.tight_layout()

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
    w_real : np.ndarray(dtype=float, ndim=1)
        Real part of the eigenvalues as a numpy array.
    w_imag : np.ndarray(dtype=float, ndim=1)
        Imaginary part of the eigenvalues as a numpy array.
    """

    def __init__(self, dataset, figsize, custom_figure, **kwargs):
        super().__init__(
            figure_type="single-spectrum", figsize=figsize, custom_figure=custom_figure
        )
        self.dataset = dataset
        self.kwargs = kwargs

        self.marker = "."
        self.color = "blue"
        self.markersize = 6
        self.alpha = 0.8

        self.w_real = self.dataset.eigenvalues.real
        self.w_imag = self.dataset.eigenvalues.imag
        self._add_spectrum()

    def _add_spectrum(self):
        (spectrum_point,) = self.ax.plot(
            self.w_real * self.x_scaling,
            self.w_imag * self.y_scaling,
            marker=self.marker,
            color=self.color,
            markersize=self.markersize,
            alpha=self.alpha,
            linestyle="None",
            **self.kwargs,
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
        c_handler : `~pylbo.continua.ContinuaHandler` instance
            The `ContinuaHandler` instance used to plot the continua.
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
            if self._c_handler.check_if_collapsed(continuum=continuum):
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
                    facecolor=color, alpha=self._c_handler.alpha_region, label=key
                )
                if key == "thermal":
                    item = self.ax.axhspan(
                        min_value * self.y_scaling, max_value * self.y_scaling, **props
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
        # this creates and sets self._ef_ax
        self._add_eigenfunction_axes()
        if self._ef_handler is None:
            self._ef_handler = EigenfunctionHandler(self.dataset, self._ef_ax)
        # connect everything
        super().add_eigenfunctions()


class MultiSpectrumPlot(SpectrumFigure):
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
        self.x_scaling = np.ones_like(self.dataseries)
        self.y_scaling = np.ones_like(self.dataseries)
        self.kwargs = kwargs
        self._add_spectrum()

    def _validate_xdata(self, xdata):
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
        ydata_values = np.array(
            [ds.eigenvalues ** self._w_pow for ds in self.dataseries]
        )
        ydata_values[np.where(np.abs(ydata_values) < 1e-12)] = np.nan
        if self.use_real_parts:
            return ydata_values.real
        else:
            return ydata_values.imag

    def set_x_scaling(self, x_scaling):
        if isinstance(x_scaling, (int, float, complex)):
            x_scaling = np.ones_like(self.dataseries) * x_scaling
        super().set_x_scaling(x_scaling)

    def set_y_scaling(self, y_scaling):
        if isinstance(y_scaling, (int, float, complex)):
            y_scaling = np.ones_like(self.dataseries) * y_scaling
        super().set_y_scaling(y_scaling)

    def _add_spectrum(self):
        """
        Draw method, creates the spectrum.
        """
        props = copy(self.kwargs)
        marker = props.pop("marker", ".")
        color = props.pop("color", "blue")
        markersize = props.pop("markersize", 6)
        alpha = props.pop("alpha", 0.8)
        for i, ds in enumerate(self.dataseries):
            (spectrum_point,) = self.ax.plot(
                self.xdata[i] * np.ones_like(self.ydata[i]) * self.x_scaling[i],
                self.ydata[i] * self.y_scaling[i],
                marker=marker,
                color=color,
                markersize=markersize,
                alpha=alpha,
                linestyle="None",
                **props,
            )
            add_pickradius_to_item(item=spectrum_point, pickradius=10)
            # set dataset associated with this line of points
            setattr(spectrum_point, "dataset", ds)
        self.ax.axhline(y=0, linestyle="dotted", color="grey", alpha=0.3)
        self.ax.axvline(x=0, linestyle="dotted", color="grey", alpha=0.3)

    def add_continua(self, interactive=True):
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
            continuum = self.dataseries.continua[key]
            min_values = (
                np.array([np.min(c_ds) for c_ds in continuum]) ** self._w_pow
                * self.y_scaling
            )
            max_values = (
                np.array([np.max(c_ds) for c_ds in continuum]) ** self._w_pow
                * self.y_scaling
            )
            # skip if continua are all zero
            if all(min_values == 0) and all(max_values == 0):
                continue
            # when continua are collapsed then min = max and we draw a line instread
            if all(abs(min_values - max_values) < 1e-12):
                (item,) = self.ax.plot(
                    self.xdata,
                    min_values,
                    color=color,
                    alpha=self._c_handler.alpha_point,
                    label=key,
                )
            else:
                item = self.ax.fill_between(
                    self.xdata,
                    min_values,
                    max_values,
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
        self._add_eigenfunction_axes()
        if self._ef_handler is None:
            self._ef_handler = EigenfunctionHandler(self.dataseries, self._ef_ax)
        # connect everything
        super().add_eigenfunctions()
