from functools import wraps
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.axes import Axes as mpl_axes
from matplotlib.figure import Figure as mpl_fig
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

    figure_stack = dict()

    def __init__(
        self,
        figure_type: str,
        figsize: tuple[int, int] = None,
        custom_figure: tuple[mpl_fig, mpl_axes] = None,
    ):
        self.fig, self.ax = self._create_axes_instances(
            custom_figure, figure_type, figsize
        )
        self.x_scaling = 1.0
        self.y_scaling = 1.0
        self._mpl_callbacks = []

    @classmethod
    def _generate_figure_id(cls, figlabel: str) -> str:
        """
        Generates a unique figure id.

        Parameters
        ----------
        figlabel : str
            The label of the figure.

        Returns
        -------
        figure_id : str
            The unique figure id of the form "figure_type-x" where x is an integer.
        """
        # count occurences of this type of id in the list
        occurences = sum(figlabel in fig_id for fig_id in list(cls.figure_stack.keys()))
        return f"{figlabel}-{1 + occurences}"

    def _create_axes_instances(
        self,
        custom_figure: tuple[mpl_fig, mpl_axes],
        figlabel: str,
        figsize: tuple[int, int],
    ) -> tuple[mpl_fig, mpl_axes]:
        """
        Creates the figure and axes instances.

        Parameters
        ----------
        custom_figure : tuple[~matplotlib.figure.Figure, ~matplotlib.axes.Axes]
            A custom figure to use, in the form (fig, ax) corresponding to the figure
            and axis objects from matplotlib.
        figlabel : str
            The label of the figure.
        figsize : tuple[int, int]
            The size of the figure, default is (10, 6).

        Returns
        -------
        fig : ~matplotlib.figure.Figure
            The figure on which to draw.
        ax : ~matplotlib.axes.Axes
            The axes on which to draw.
        """
        if custom_figure is not None:
            fig, ax = custom_figure
            figsize = tuple(fig.get_size_inches())
            fig_id = fig.get_label()
        else:
            figsize = figsize if figsize is not None else (10, 6)
            fig_id = self._generate_figure_id(figlabel)
            fig = plt.figure(fig_id, figsize=figsize)
            ax = fig.add_subplot(111)
        self.figsize = figsize
        self.figure_id = fig_id
        self.__class__.figure_stack.update({fig_id: self})
        return fig, ax

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

        _geometry = transform_to_numpy(get_axis_geometry(ax))
        subplot_index = _geometry[-1]
        gspec_outer = gridspec.GridSpec(*_geometry[0:2], figure=self.fig)
        gspec_inner = gridspec.GridSpecFromSubplotSpec(
            *subplot_geometry, subplot_spec=gspec_outer[subplot_index]
        )
        ax.set_subplotspec(gspec_inner[old_new_position[0]])
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
        """Shows the current figure"""
        plt.figure(self.figure_id)

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
