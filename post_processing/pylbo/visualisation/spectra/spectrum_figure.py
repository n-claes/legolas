from copy import copy

from pylbo.visualisation.figure_manager import FigureWindow


class SpectrumFigure(FigureWindow):
    """Superclass of both single and multispectra."""

    def __init__(self, figure_type, figsize, custom_figure):
        super().__init__(
            figure_type=figure_type, figsize=figsize, custom_figure=custom_figure
        )
        self._c_handler = None
        self._ef_handler = None
        self._ef_ax = None
        self._def_handler = None
        self._def_ax = None

        self.plot_props = None
        self.marker = None
        self.color = None
        self.markersize = None
        self.alpha = None

    def draw(self):
        """
        Called on plot refreshing, the super call clears the figure and then redraws
        everything.
        """
        super().draw()
        self._add_spectrum()
        if self._c_handler is not None:
            self.add_continua(self._c_handler.interactive)
        if self._ef_handler is not None:
            self.add_eigenfunctions()
        if self._def_handler is not None:
            self.add_derived_eigenfunctions()

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
        Adds eigenfunctions to the spectrum, overridden in subclasses.
        """
        pass

    def add_derived_eigenfunctions(self):
        """
        Adds derived eigenfunctions to the spectrum, overridden in subclasses.
        """
        pass

    @property
    def c_handler(self):
        """Property, returns the continua handler."""
        return self._c_handler

    @property
    def ef_handler(self):
        """Property, returns the eigenfunction handler."""
        return self._ef_handler

    @property
    def ef_ax(self):
        """Property, returns the eigenfunction axes."""
        return self._ef_ax
