import numpy as np
from matplotlib.collections import PathCollection
from pylbo.utilities.logger import pylboLogger
from pylbo.utilities.toolbox import add_pickradius_to_item
from pylbo.data_containers import LegolasDataSet, LegolasDataSeries
from pylbo.exceptions import PostprocessedNotPresent


class PostprocessedHandler:
    """
    Main handler for post-processed quantities.
    """

    def __init__(self, data, pp_ax):
        self.data = data
        self.pp_ax = pp_ax
        # check presence of post-processed quantities
        error_msg = "None of the given datfiles has post-processed quantities written to it."
        if isinstance(self.data, LegolasDataSeries):
            if not any(self.data.pp_written):
                raise PostprocessedNotPresent(error_msg)
        else:
            if not self.data.pp_written:
                raise PostprocessedNotPresent(error_msg)

        # holds the points that are currently selected in the form of a "double" dict:
        # {"ds instance" : {"index" : line2D instance}}
        self._selected_idxs = {}
        self._use_real_part = True
        # index for the selected quantity, cfr. state vector (so 0 = "rho")
        self._selected_pp_name_idx = 0
        self._pp_names = self.data.pp_names
        # flag to check if points are marked
        self._no_pp_is_transparent = False
        self._unmarked_alpha = None

    def on_point_pick(self, event):
        """
        Determines what happens when a pickable artist is selected.

        Parameters
        ----------
        event : ~matplotlib.backend_bases.PickEvent
            The pick event.
        """
        # comparing spectra: prevents double event triggering
        if hasattr(self, "associated_ds_ax"):
            if not event.mouseevent.inaxes == self.associated_ds_ax:
                return
        artist = event.artist
        # if artist is a legend item, return (this attribute has been set manually)
        if hasattr(artist, "is_legend_item"):
            return
        # retrieve figure and axis for this artist, save limits to prevent zoom reset
        fig, ax = artist.figure, artist.axes
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        # Retrieves the dataset associated with the current points.
        # This attribute has been set when the spectrum was added to the plot.
        if not hasattr(artist, "dataset"):
            return
        associated_ds = artist.dataset
        # This retrieves the indices of the clicked points. Multiple indices are
        # possible depending on an overlapping pickradius. Look which point corresponds
        # to the smallest distance to the mouse click.
        idxs = event.ind
        if isinstance(artist, PathCollection):
            # this is done for points drawn using scatter instead of plot
            xdata, ydata = np.split(artist.get_offsets(), [-1], axis=1)
            xdata = xdata.squeeze()
            ydata = ydata.squeeze()
        else:
            xdata = artist.get_xdata()
            ydata = artist.get_ydata()
        if len(idxs) == 1:
            idx = idxs[0]
        else:
            mouse_x = event.mouseevent.xdata
            mouse_y = event.mouseevent.ydata
            distances = (mouse_x - xdata[idxs]) ** 2 + (mouse_y - ydata[idxs]) ** 2
            idx = idxs[distances.argmin()]
        xdata = xdata[idx]
        ydata = ydata[idx]

        # if the selected point has no post-processed quantities, do nothing
        if not associated_ds.pp_written:
            pylboLogger.warn(
                f"the selected point at ({xdata: 3.2e}, {ydata: 3.2e}) "
                f"has no post-processed quantities associated with it. "
                f"(corresponding dataset: {associated_ds.datfile.stem})"
            )
            return

        # handle left clicking
        if event.mouseevent.button == 1:
            # skip if point index is already in list
            if str(idx) in self._selected_idxs.get(associated_ds, {}).keys():
                return
            (marked_point,) = ax.plot(
                xdata,
                ydata,
                "x",
                markersize=8,
                label="marked_point",
                alpha=0.7,
                markeredgewidth=3,
            )
            add_pickradius_to_item(item=marked_point, pickradius=1)
            # get items corresponding to this ds
            items = self._selected_idxs.get(associated_ds, {})
            items.update({f"{idx}": marked_point})
            self._selected_idxs.update({associated_ds: items})

        # handle right clicking
        elif event.mouseevent.button == 3:
            # remove selected index from list
            selected_artist = self._selected_idxs.get(associated_ds, {}).pop(
                str(idx), None
            )
            if selected_artist is not None:
                selected_artist.remove()
                # if no items remaining for this ds, remove key
                if len(self._selected_idxs[associated_ds]) == 0:
                    self._selected_idxs.pop(associated_ds)

        # this fixes a zoom reset when picking the figure
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        fig.canvas.draw()

    def on_key_press(self, event):
        """
        Determines what happens when a key is pressed.

        Parameters
        ----------
        event : ~matplotlib.backend_bases.KeyEvent
            The key event.
        """
        # pressing "m" marks points without post-processed quantities
        if event.key == "m":
            self.mark_points_without_ppqs()
        # pressing "d" clears figure and selection
        if event.key == "d":
            # retrieve (index, item) dicts for each dataset
            for ds_artist_dict in self._selected_idxs.values():
                for artist in ds_artist_dict.values():
                    artist.remove()
            self._selected_idxs.clear()
            self.pp_ax.clear()
        # pressing "i" switches between real and imaginary parts
        if event.key == "i":
            self._use_real_part = not self._use_real_part
        # pressing "up" or "down" cycles through the post-processed quantities
        if event.key == "up":
            self._select_next_ppq()
        if event.key == "down":
            self._select_previous_ppq()
        # pressing one of these keys updates the plot
        if event.key in ("enter", "up", "down", "i", "t"):
            self.update_plot()
        event.canvas.draw()

    def update_plot(self):
        """
        Updates the plot when an event is triggered, clears and then redraws
        the post-processed quantities. Rescaling of the axes is done automatically.
        """
        # do nothing if nothing is selected
        if not self._selected_idxs:
            self.pp_ax.clear()
            return
        self.pp_ax.clear()
        pp_name = self._pp_names[self._selected_pp_name_idx]
        for ds, idxs_dict in self._selected_idxs.items():
            idxs = np.array([int(idx) for idx in idxs_dict.keys()])
            pp_container = ds.get_postprocessed(ev_idxs=idxs)

            for ev_idx, ppqs in zip(idxs, pp_container):
                pp = ppqs.get(pp_name)
                if self._use_real_part:
                    pp = pp.real
                else:
                    pp = pp.imag
                # create label associated with ds if it's a series
                label = rf"$\omega_{{{ev_idx}}}$ = {ppqs.get('eigenvalue'):2.3e}"
                if isinstance(self.data, LegolasDataSeries):
                    label = "".join([f"{ds.datfile.stem} | ", label])
                # get color of selected point
                color = self._selected_idxs.get(ds).get(str(ev_idx)).get_color()
                self.pp_ax.plot(ds.ef_grid, pp, color=color, label=label)
        self.pp_ax.axhline(y=0, linestyle="dotted", color="grey")
        if isinstance(self.data, LegolasDataSet):
            self.pp_ax.axvline(x=self.data.x_start, linestyle="dotted", color="grey")
        self.pp_ax.set_title(self._get_title(pp_name))
        self.pp_ax.legend(loc="best")
        # updates gridspec
        self.pp_ax.figure.tight_layout()

    def _select_next_ppq(self):
        """Increments the post-processed quantity index by one."""
        self._selected_pp_name_idx += 1
        if self._selected_pp_name_idx > len(self._pp_names) - 1:
            self._selected_pp_name_idx = 0

    def _select_previous_ppq(self):
        """Decrements the post-processed quantity index by one."""
        self._selected_pp_name_idx -= 1
        if self._selected_pp_name_idx < 0:
            self._selected_pp_name_idx = len(self._pp_names) - 1

    def _get_title(self, pp_name):
        """
        Creates the title of the post-processed quantity plot.
        LaTeX formatting is used and numbers are replaced with the
        corresponding suffix: :math:`(1, 2, 3)`
        becomes :math:`(x, y, z)` or :math:`(r, \\theta, z)`.

        Parameters
        ----------
        pp_name : str
            The name of the post-processed quantity as present in the datfile.

        Returns
        -------
        name : str
            The 'new' name for the post-processed quantity, used as title.
        """
        part = "Re"
        if not self._use_real_part:
            part = "Im"
        suffix = ("_x", "_y", "_z")
        if self.data.geometry == "cylindrical":
            suffix = ("_r", r"_\theta", "_z")
        for i, idx in enumerate("123"):
            pp_name = pp_name.replace(idx, suffix[i])
        pp_name = pp_name.replace("div", "\\nabla\\cdot")
        pp_name = pp_name.replace("curl", "\\nabla\\times")
        pp_name = pp_name.replace("para", "\\parallel")
        pp_name = pp_name.replace("perp", "\\perp")
        name = fr"{part}(${pp_name}$)"
        return name

    def mark_points_without_ppqs(self):
        """
        For dataseries, it is possible that not all datasets in the series
        have post-processed quantities associated with them. This routine will
        toggle a change in the opacity value for datapoints with no post-processed
        quantities, so they are clearly distinguishable from those who do have them.
        """
        if not isinstance(self.data, LegolasDataSeries):
            return
        if all(self.data.pp_written):
            return
        self._no_pp_is_transparent = not self._no_pp_is_transparent
        for ax in self.pp_ax.figure.get_axes():
            for child in ax.get_children():
                # the ones with this attribute are all the vertical rows of datapoints
                if hasattr(child, "dataset"):
                    if self._unmarked_alpha is None:
                        self._unmarked_alpha = child.get_alpha()
                    if not child.dataset.pp_written:
                        if self._no_pp_is_transparent:
                            child.set_alpha(0)
                        else:
                            child.set_alpha(self._unmarked_alpha)
