import numpy as np
from pylbo.utilities.logger import pylboLogger
from pylbo.utilities.toolbox import add_pickradius_to_item
from pylbo.data_containers import LegolasDataSet, LegolasDataSeries
from pylbo.exceptions import EigenfunctionsNotPresent


class EigenfunctionHandler:
    """
    Main handler for eigenfunctions.
    """

    def __init__(self, data, ef_ax):
        self.data = data
        self.ef_ax = ef_ax
        # check presence of eigenfunctions
        error_msg = "None of the given datfiles has eigenfunctions written to it."
        if isinstance(self.data, LegolasDataSeries):
            if not any(self.data.efs_written):
                raise EigenfunctionsNotPresent(error_msg)
        else:
            if not self.data.efs_written:
                raise EigenfunctionsNotPresent(error_msg)

        # holds the points that are currently selected in the form of a "double" dict:
        # {"ds instance" : {"index" : line2D instance}}
        self._selected_idxs = {}
        self._use_real_part = True
        # index for the selected eigenfunction, cfr. state vector (so 0 = "rho")
        self._selected_ef_name_idx = 0
        self._ef_names = self.data.ef_names
        # if True, retransforms the eigenfunctions to e.g. rv_r in cylindrical geometry
        self._retransform_efs = False
        # flag to check if points are marked
        self._no_efs_is_transparent = False
        self._unmarked_alpha = None

    def on_point_pick(self, event):
        """
        Determines what happens when a pickable artist is selected.

        Parameters
        ----------
        event : ~matplotlib.backend_bases.PickEvent
            The pick event.
        """
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
        associated_ds = artist.dataset
        # This retrieves the indices of the clicked points. Multiple indices are
        # possible depending on an overlapping pickradius. Look which point corresponds
        # to the smallest distance to the mouse click.
        idxs = event.ind
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

        # if the selected point has no eigenfunctions, do nothing
        if not associated_ds.efs_written:
            pylboLogger.warn(
                f"the selected point at ({xdata: 3.2e}, {ydata: 3.2e}) "
                f"has no eigenfunctions associated with it. "
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
        # pressing "m" marks points without eigenfunctions
        if event.key == "m":
            self.mark_points_without_eigenfunctions()
        # pressing "d" clears figure and selection
        if event.key == "d":
            # retrieve (index, item) dicts for each dataset
            for ds_artist_dict in self._selected_idxs.values():
                for artist in ds_artist_dict.values():
                    artist.remove()
            self._selected_idxs.clear()
            self.ef_ax.clear()
        # pressing "i" switches between real and imaginary parts
        if event.key == "i":
            self._use_real_part = not self._use_real_part
        # pressing "up" or "down" cycles through the eigenfunctions
        if event.key == "up":
            self._select_next_eigenfunction()
        if event.key == "down":
            self._select_previous_eigenfunction()
        # pressing "t" retransforms the eigenfunctions
        if event.key == "t":
            self._retransform_efs = not self._retransform_efs
        # pressing one of these keys updates the plot
        if event.key in ("enter", "up", "down", "i", "t"):
            self.update_plot()
        event.canvas.draw()

    def update_plot(self):
        """
        Updates the plot when an event is triggered, clears and then redraws
        the eigenfunctions. Rescaling of the axes is done automatically.
        """
        # do nothing if nothing is selected
        if not self._selected_idxs:
            self.ef_ax.clear()
            return
        self.ef_ax.clear()
        ef_name = self._ef_names[self._selected_ef_name_idx]
        for ds, idxs_dict in self._selected_idxs.items():
            idxs = np.array([int(idx) for idx in idxs_dict.keys()])
            ef_container = ds.get_eigenfunctions(ev_idxs=idxs)

            for ev_idx, efs in zip(idxs, ef_container):
                ef = efs.get(ef_name)
                if self._use_real_part:
                    ef = ef.real
                else:
                    ef = ef.imag
                # check for retransform
                if (
                    self._retransform_efs
                    and ef_name in ("rho", "v1", "v3", "T", "a2")
                    and self.data.geometry == "cylindrical"
                ):
                    ef = ef * ds.ef_grid
                # create label associated with ds if it's a series
                label = rf"$\omega_{{{ev_idx}}}$ = {efs.get('eigenvalue'):2.3e}"
                if isinstance(self.data, LegolasDataSeries):
                    label = "".join([f"{ds.datfile.stem} | ", label])
                # get color of selected point
                color = self._selected_idxs.get(ds).get(str(ev_idx)).get_color()
                self.ef_ax.plot(ds.ef_grid, ef, color=color, label=label)
        self.ef_ax.axhline(y=0, linestyle="dotted", color="grey")
        if isinstance(self.data, LegolasDataSet):
            self.ef_ax.axvline(x=self.data.x_start, linestyle="dotted", color="grey")
        self.ef_ax.set_title(self._get_title(ef_name))
        self.ef_ax.legend(loc="best")
        # updates gridspec
        self.ef_ax.figure.tight_layout()

    def _select_next_eigenfunction(self):
        """Increments the eigenfunction index by one."""
        self._selected_ef_name_idx += 1
        if self._selected_ef_name_idx > len(self._ef_names) - 1:
            self._selected_ef_name_idx = 0

    def _select_previous_eigenfunction(self):
        """Decrements the eigenfunction index by one."""
        self._selected_ef_name_idx -= 1
        if self._selected_ef_name_idx < 0:
            self._selected_ef_name_idx = len(self._ef_names) - 1

    def _get_title(self, ef_name):
        """
        Creates the title of the eigenfunction plot.
        If the eigenfunction is retransformed an 'r' is prepended if appropriate,
        along with Re/Im depending on the real/imaginary part shown.
        Additionally, LaTeX formatting is used and numbers are replaced with the
        corresponding suffix: :math:`(1, 2, 3)`
        becomes :math:`(x, y, z)` or :math:`(r, \\theta, z)`.

        Parameters
        ----------
        ef_name : str
            The name of the eigenfunction as present in the datfile.

        Returns
        -------
        name : str
            The 'new' name for the eigenfunction, used as title.
        """
        part = "Re"
        if not self._use_real_part:
            part = "Im"
        suffix = ("_x", "_y", "_z")
        if self.data.geometry == "cylindrical":
            suffix = ("_r", r"_\theta", "_z")
            if ef_name in ("rho", "v1", "v3", "T", "a2") and self._retransform_efs:
                ef_name = f"r{ef_name}"
        for i, idx in enumerate("123"):
            ef_name = ef_name.replace(idx, suffix[i])
        ef_name = ef_name.replace("rho", r"\rho")
        name = f"{part}(${ef_name}$) eigenfunction"
        return name

    def mark_points_without_eigenfunctions(self):
        """
        For dataseries, it is possible that not all datasets in the series
        have eigenfunctions associated with them. This routine will toggle a change
        in the opacity value for datapoints with no eigenfunctions, so they are
        clearly distinguishable from those who do have them.
        """
        if not isinstance(self.data, LegolasDataSeries):
            return
        if all(self.data.efs_written):
            return
        self._no_efs_is_transparent = not self._no_efs_is_transparent
        for ax in self.ef_ax.figure.get_axes():
            for child in ax.get_children():
                # the ones with this attribute are all the vertical rows of datapoints
                if hasattr(child, "dataset"):
                    if self._unmarked_alpha is None:
                        self._unmarked_alpha = child.get_alpha()
                    if not child.dataset.efs_written:
                        if self._no_efs_is_transparent:
                            child.set_alpha(0)
                        else:
                            child.set_alpha(self._unmarked_alpha)
