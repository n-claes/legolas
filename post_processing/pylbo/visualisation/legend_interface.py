import matplotlib.collections as mpl_collections
import matplotlib.lines as mpl_lines
import numpy as np
from pylbo.utilities.toolbox import add_pickradius_to_item


class LegendHandler:
    def __init__(self, interactive):
        self.legend = None
        self.alpha_point = 0.8
        self.alpha_region = 0.2
        self.alpha_hidden = 0.05

        self.marker = "p"
        self.markersize = 64
        self.pickradius = 10
        self.linewidth = 2
        self.legend_properties = {}

        self.interactive = interactive
        self._drawn_items = []
        self._legend_mapping = {}

    def on_legend_pick(self, event):
        """
        Determines what happens when the legend gets picked.

        Parameters
        ----------
        event : `~matplotlib.backend_bases.PickEvent`
            The matplotlib pick event.
        """
        artist = event.artist
        if artist not in self._legend_mapping:
            return
        drawn_item = self._legend_mapping.get(artist)
        visible = not drawn_item.get_visible()
        drawn_item.set_visible(visible)
        if visible:
            if isinstance(artist, (mpl_collections.PathCollection, mpl_lines.Line2D)):
                artist.set_alpha(self.alpha_point)
            else:
                artist.set_alpha(self.alpha_region)
        else:
            artist.set_alpha(self.alpha_hidden)
        artist.figure.canvas.draw()

    def make_legend_pickable(self):
        """
        Makes the legend pickable, only used if interactive.
        """
        legend_handles = self.legend.legendHandles
        handle_labels = [handle.get_label() for handle in legend_handles]
        # we need a mapping of the legend item to the actual item that was drawn
        for i, drawn_item in enumerate(self._drawn_items):
            # TODO: for some reason fill_between returns empty handles such that
            #       the code below errors out. The try-except clause is an attempt
            #       to fix it and _should_ work. This relies on the ordening of the
            #       legend items when they are drawn, but should be thorougly tested.
            try:
                idx = handle_labels.index(drawn_item.get_label())
                legend_item = self.legend.legendHandles[idx]
            except ValueError:
                idx = i
                legend_item = self.legend.legendHandles[idx]
                # fix empty label
                legend_item.set_label(drawn_item.get_label())
            if not drawn_item.get_label() == legend_item.get_label():
                raise ValueError(
                    f"something went wrong in mapping legend items to drawn items. \n"
                    f"Tried to map {legend_item} (label '{legend_item.get_label()}')"
                    f" to {drawn_item} (label '{drawn_item.get_label()}') \n"
                )
            add_pickradius_to_item(item=legend_item, pickradius=self.pickradius)
            legend_item.set_alpha(self.alpha_hidden)
            # add an attribute to this artist to tell it's from a legend
            setattr(legend_item, "is_legend_item", True)
            # we make the regions invisible until clicked
            drawn_item.set_visible(False)
            self._legend_mapping[legend_item] = drawn_item

    def add(self, item):
        """
        Adds an item to the list of drawn items on the canvas.

        Parameters
        ----------
        item : object
            A single object, usually a return from the matplotlib plot or scatter
            methods.
        """
        if isinstance(item, (list, np.ndarray, tuple)):
            raise ValueError("object expected, not something list-like")
        self._drawn_items.append(item)
