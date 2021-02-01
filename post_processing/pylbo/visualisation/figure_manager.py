import matplotlib.pyplot as plt
from pathlib import Path
from functools import wraps
from pylbo.utilities.logger import pylboLogger


def refresh_plot(f):
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
        when calling plt.show(). If `False`, the figures in the dictionary will remain
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
        figure : `~pylbo.visualisation.figure_manager.PlotWindow` instance
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
        figure : `~pylbo.visualisation.figure_manager.PlotWindow` instance
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

    def disable_stack(self):
        """
        This disables the figurestack by closing the figures: this prevents them
        from showing when calling plt.show(). Later on we reconstruct the
        figuremanagers.
        """
        [figure.disable() for figure in self.figure_list]
        self.stack_is_enabled = False

    def enable_stack(self):
        """
        This enables the figurestack so figures are again accessible through
        the plt.show() methods. When we lock the stack those figures are closed,
        which essentially destroys the figuremanager but keeps the figure references
        alive. Here we simply reconstruct the figuremanager for every disabled figure
        in the stack, allowing for multiple calls to show(), save(), etc.
        """
        [figure.enable() for figure in self.figure_list]
        self.stack_is_enabled = True

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
    fig : matplotlib.figure.Figure
        The figure on which to draw.
    ax : matplotlib.axes.Axes
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
    enabled : bool
        If `False` (default), the current figure will not be plotted at plt.show().
    """

    figure_stack = FigureContainer()

    def __init__(self, figure_type, figsize=None, custom_figure=None):
        if custom_figure is not None:
            plt.close(custom_figure[0])
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
        self.x_scaling = 1
        self.y_scaling = 1
        self.enabled = False
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

    def disable(self):
        """Disables the current figure."""
        self.enabled = False
        plt.close(self.figure_id)

    def enable(self):
        """Enables the current figure, reconstructs the figuremanagers"""
        self.enabled = True
        manager = plt.figure(self.figure_id, self.figsize).canvas.manager
        manager.canvas.figure = self.fig
        self.fig.set_canvas(manager.canvas)

    def connect_callbacks(self):
        for callback in self._mpl_callbacks:
            self.fig.canvas.mpl_connect(callback["kind"], callback["method"])

    def disconnect_callbacks(self):
        for callback in self._mpl_callbacks:
            self.fig.canvas.mpl_disconnect(callback["cid"])

    def draw(self):
        self.ax.cla()
        self.disconnect_callbacks()

    @refresh_plot
    def set_x_scaling(self, x_scaling):
        self.x_scaling = x_scaling

    @refresh_plot
    def set_y_scaling(self, y_scaling):
        self.y_scaling = y_scaling

    def show(self):
        """Shows the current figure."""
        self.__class__.figure_stack.disable_stack()
        self.enable()
        self.connect_callbacks()
        plt.show()

    @classmethod
    def showall(cls):
        """Shows all active figures at once through a call to plt.show()."""
        cls.figure_stack.disable_stack()
        for figure in cls.figure_stack.figure_list:
            figure.enable()
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
            Default keyword arguments passed to `matplotlib.pyplot.savefig`.
        """
        filepath = Path(filename).resolve()
        self.fig.savefig(filepath, **kwargs)
        pylboLogger.info(f"figure saved to {filepath}")
