import matplotlib.pyplot as plt


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
        if self.stack_is_enabled:
            return
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
        if not self.stack_is_enabled:
            return
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
        if not self.enabled:
            return
        self.enabled = False
        plt.close(self.figure_id)

    def enable(self):
        """Enables the current figure, reconstructs the figuremanagers"""
        if self.enabled:
            return
        self.enabled = False
        manager = plt.figure(self.figure_id, self.figsize).canvas.manager
        manager.canvas.figure = self.fig
        self.fig.set_canvas(manager.canvas)

    def draw(self):
        """Method to draw, should be overriden by subclass."""
        pass

    def add_continua(self):
        """Method to add continua, should be overridden by subclass."""
        pass

    def add_eigenfunctions(self):
        """Method to add eigenfunctions, should be overridden by subclass."""
        pass

    def show(self):
        """Shows the current figure."""
        self.__class__.figure_stack.disable_stack()
        self.enable()
        self.draw()
        plt.show()
        self.__class__.figure_stack.enable_stack()

    @classmethod
    def showall(cls, reconstruct=False):
        """
        Shows all active figures at once through a call to plt.show().
        Unless `reconstruct=True` this is final: all figures and figuremanagers
        will be destroyed and will not be reconstructed.

        Parameters
        ----------
        reconstruct : bool
            If `True`, reconstructs all figuremanagers, allowing for multiple
            calls to the show routines.

        Raises
        ------
        ValueError
            If no active figures are present.
        ValueError
            If the figure stack is disabled.
        """
        if cls.figure_stack.is_empty:
            raise ValueError("No active figures present, stack is empty.")
        if not cls.figure_stack.stack_is_enabled:
            raise ValueError("Figure stack is disabled.")
        for figure in cls.figure_stack.figure_list:
            figure.draw()
        plt.show()
        if reconstruct:
            cls.figure_stack.enable_stack()
        else:
            cls.figure_stack.clear()

    def save(self, filename):
        """
        Saves the current figure.

        Parameters
        ----------
        filename : str, ~os.PathLike
            The filename to which the current figure is saved.

        Raises
        ------
        ValueError
            If the current figure is disabled.
        """
        if not self.enabled:
            raise ValueError("This figure is disabled.")
        self.__class__.figure_stack.disable_stack()
        self.enable()
        self.fig.save(filename)
        self.__class__.figure_stack.enable_stack()

    def clear(self):
        """Clears the current axes."""
        self.ax.clear()
