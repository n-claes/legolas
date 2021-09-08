import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from functools import wraps
from pylbo.utilities.logger import pylboLogger
from pylbo.utilities.toolbox import get_axis_geometry, transform_to_numpy


def refresh_plot(f):
    """
    Simple decorator, when a routine is wrapped with this the plot will be
    cleared and redrawn, useful for when the scaling is changed or artists are
    added/removed.
    """

    @wraps(f)
    def refresh(*args, **kwargs):
        f(*args, **kwargs)
        window = args[0]
        window.draw()
        return f

    return refresh


class FigureContainer(dict):
    """
    A special dictionary containing the currently active figures.

    Attributes
    ----------
    stack_is_enabled : bool
        If `True` (default), the dictionary is unlocked and figures are drawn
        when calling `plt.show()`.
        If `False`, the figures in the dictionary will remain
        closed and will not be drawn.
    """

    def __init__(self):
        super().__init__()
        self.stack_is_enabled = True

    def add(self, figure):
        """
        Adds a new figure to the stack.

        Parameters
        ----------
        figure : ~pylbo.visualisation.figure_manager.FigureWindow
            The figure to add.

        Raises
        -------
        ValueError
            If the figure id is already present in the stack.
        """
        if figure.figure_id in self.figure_id_list:
            raise ValueError(
                f"id = '{figure.figure_id}' already in existing list: "
                f"{self.figure_id_list}"
            )
        self.update({figure.figure_id: figure})

    def pop(self, figure_id):
        """
        Removes and returns the figure corresponding to the given id from the stack.

        Parameters
        ----------
        figure_id : str
            The figure id, corresponds to the dictionary key.

        Returns
        -------
        figure : ~pylbo.visualisation.figure_manager.FigureWindow
            The figure corresponding to `figure_id`.
        """
        self._validate_figure_id(figure_id)
        plt.close(figure_id)
        return super().pop(figure_id)

    def _validate_figure_id(self, figure_id):
        """
        Checks if the current id is present, needed to avoid matplotlib errors
        when trying to close.

        Parameters
        ----------
        figure_id : str
            The figure id, corresponds to the dictionary key.

        Raises
        -------
        ValueError
            If the given id is not present in the list.
        """
        if figure_id not in self.figure_id_list:
            raise ValueError(
                f"id='{figure_id}' not in existing list {self.figure_id_list}"
            )

    @property
    def number_of_figures(self):
        """Returns the total number of figures in the stack."""
        return len(self)

    @property
    def figure_id_list(self):
        """Returns the list of figure ids in the stack."""
        return list(self.keys())

    @property
    def figure_list(self):
        """Returns the list of figures in the stack."""
        return list(self.values())

    @property
    def is_empty(self):
        """Returns `True` if there are no active figures."""
        return len(self) == 0


class FigureWindow:
    """
    Class to handle drawing different types of plots.

    Parameters
    ----------
    figure_type : str
        The type of figure to create, links to the figure id.
    figsize : tuple
        The figuresize as a matplotlib (x, x) tuple.
    custom_figure : tuple
        A custom figure to use, in the form (fig, ax) corresponding to the figure
        and axis objects from matplotlib.

    Attributes
    ----------
    fig : ~matplotlib.figure.Figure
        The figure on which to draw.
    ax : ~matplotlib.axes.Axes
        The axes on which to draw.
    custom_figure : tuple
        The (fig, ax) custom figure that was passed, `None` by default.
    figsize : tuple
        The size of the figure.
    figure_id : str
        A unique id associated with the figure.
    x_scaling : int, float, complex, np.ndarray
        Scaling to apply to the x-axis.
    y_scaling : int, float, complex, np.ndarray
        Scaling to apply to the y-axis.
    """

    figure_stack = FigureContainer()

    def __init__(self, figure_type, figsize=None, custom_figure=None):
        if custom_figure is not None:
            self.fig, self.ax = custom_figure
            self.figsize = tuple(self.fig.get_size_inches())
            if self.fig.get_label() == "":
                self.figure_id = self._generate_figure_id(figure_type)
            else:
                self.figure_id = self._generate_figure_id(self.fig.get_label())
        else:
            self.figsize = figsize
            if self.figsize is None:
                self.figsize = (12, 8)
            self.figure_id = self._generate_figure_id(figure_type)
            self.fig = plt.figure(self.figure_id, figsize=self.figsize)
            self.ax = self.fig.add_subplot(111)
        self.x_scaling = 1.0
        self.y_scaling = 1.0
        self._mpl_callbacks = []
        self.__class__.figure_stack.add(self)

    @classmethod
    def _generate_figure_id(cls, figure_type):
        """
        Generates a unique figure id.

        Parameters
        ----------
        figure_type : str
            The type of figure to create

        Returns
        -------
        figure_id : str
            The unique figure id of the form "figure_type-x" where x is an integer.
        """
        # count occurences of this type of id in the list
        occurences = sum(
            figure_type in fig_id for fig_id in cls.figure_stack.figure_id_list
        )
        suffix = 1 + occurences
        return f"{figure_type}-{suffix}"

    def connect_callbacks(self):
        """Connects all callbacks to the canvas"""
        for callback in self._mpl_callbacks:
            self.fig.canvas.mpl_connect(callback["kind"], callback["method"])

    def disconnect_callbacks(self):
        """Disconnects all callbacks from the canvas"""
        for callback in self._mpl_callbacks:
            self.fig.canvas.mpl_disconnect(callback["cid"])

    def _enable_interactive_legend(self, handle):
        """
        Makes the current legendhandler interactive.

        Parameters
        ----------
        handle : ~pylbo.visualisation.legend_interface.LegendHandler
            The legendhandler (or a subclass thereof)
        """
        handle.make_legend_pickable()
        callback_kind = "pick_event"
        callback_method = handle.on_legend_pick
        callback_id = self.fig.canvas.mpl_connect(callback_kind, callback_method)
        self._mpl_callbacks.append(
            {
                "cid": callback_id,
                "kind": callback_kind,
                "method": callback_method,
            }
        )

    def _enable_interface(self, handle):
        """
        Enables the EigenfunctionInterface based on a given handler.

        Parameters
        ----------
        handle : ~pylbo.eigenfunction_interface.EigenfunctionInterface
            The handler to add to the figure.
        """
        callback_kinds = ("pick_event", "key_press_event")
        callback_methods = (handle.on_point_pick, handle.on_key_press)
        for callback_kind, callback_method in zip(callback_kinds, callback_methods):
            callback_id = self.fig.canvas.mpl_connect(callback_kind, callback_method)
            self._mpl_callbacks.append(
                {"cid": callback_id, "kind": callback_kind, "method": callback_method}
            )

    def _add_subplot_axes(self, ax, loc="right", share=None, apply_tight_layout=True):
        """
        Adds a new subplot to a given matplotlib subplot, essentially "splitting" the
        axis into two. Position and placement depend on the loc argument.
        When called on a more complex subplot layout the overall gridspec remains
        untouched, only the `ax` object has its gridspec modified.
        On return, `tight_layout()` is called by default to prevent overlapping labels.

        Parameters
        ----------
        ax : ~matplotlib.axes.Axes
            The axes object, this will be "split" and a new axes will be added
            to the figure.
        loc : str
            The location of the new subplot to add, should be "right", "left", "top"
            or "bottom". Equal to "right" by default, so the original figure is shifted
            to the left and the new axis is added on the right.
        share : str
            Can be "x", "y" or "all". This locks axes zooming between both subplots.
        apply_tight_layout : bool
            If `True` applies `tight_layout()` to the figure on return.

        Raises
        ------
        ValueError
            If the loc argument is invalid.

        Returns
        -------
        ~matplotlib.axes.Axes
            The axes instance that was added.
        """
        sharex = None
        sharey = None
        if share in ("x", "all"):
            sharex = ax
        if share in ("y", "all"):
            sharey = ax

        if loc == "right":
            subplot_geometry = (1, 2)
            old_new_position = (0, 1)
        elif loc == "left":
            subplot_geometry = (1, 2)
            old_new_position = (1, 0)
        elif loc == "top":
            subplot_geometry = (2, 1)
            old_new_position = (1, 0)
        elif loc == "bottom":
            subplot_geometry = (2, 1)
            old_new_position = (0, 1)
        else:
            raise ValueError(
                f"invalid loc={loc}, expected ['top', 'right', 'bottom', 'left']"
            )

        _geometry = transform_to_numpy(get_axis_geometry(self.ax))
        subplot_index = _geometry[-1]
        gspec_outer = gridspec.GridSpec(*_geometry[0:2], figure=self.fig)
        gspec_inner = gridspec.GridSpecFromSubplotSpec(
            *subplot_geometry, subplot_spec=gspec_outer[subplot_index]
        )
        self.ax.set_subplotspec(gspec_inner[old_new_position[0]])
        new_axis = self.fig.add_subplot(
            gspec_inner[old_new_position[1]], sharex=sharex, sharey=sharey
        )
        if apply_tight_layout:
            self.fig.tight_layout()
        return new_axis

    def draw(self):
        """
        Draws everything on the figure. This method is meant to be overriden by
        a subclass, where a `super()` call can be done to this method to clear the
        figure first.
        """
        self.ax.cla()
        self.disconnect_callbacks()

    @refresh_plot
    def set_x_scaling(self, x_scaling):
        """
        Sets the x scaling, should be overridden by a subclass.

        Parameters
        ----------
        x_scaling : int, float, complex, numpy.ndarray
            The scaling to apply to the x-axis.
        """
        self.x_scaling = x_scaling

    @refresh_plot
    def set_y_scaling(self, y_scaling):
        """
        Sets the y scaling, should be overridden by a subclass.

        Parameters
        ----------
        y_scaling : int, float, complex, numpy.ndarray
            The scaling to apply to the y-axis
        """
        self.y_scaling = y_scaling

    def show(self):
        """Shows the figures, wrapper to `showall()` for backwards compatibility."""
        self.showall()

    @classmethod
    def showall(cls):
        """Shows all active figures at once through a call to plt.show()."""
        for figure in cls.figure_stack.figure_list:
            figure.connect_callbacks()
        plt.show()

    def save(self, filename, **kwargs):
        """
        Saves the current figure.

        Parameters
        ----------
        filename : str, ~os.PathLike
            The filename to which the current figure is saved.
        kwargs
            Default keyword arguments passed to :meth:`~matplotlib.pyplot.savefig`.
        """
        filepath = Path(filename).resolve()
        self.fig.savefig(filepath, **kwargs)
        pylboLogger.info(f"figure saved to {filepath}")
