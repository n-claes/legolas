import numpy as np
from pylbo.utilities.toolbox import add_pickradius_to_item
from pylbo.utilities.datfile_utils import read_eigenfunction
from pylbo.data_containers import LegolasDataSet, LegolasDataSeries
from pylbo.exceptions import EigenfunctionsNotPresent


def get_eigenfunctions(dataset, ef_idxs):
    if not isinstance(dataset, LegolasDataSet):
        raise ValueError("get_eigenfunctions should be called on a single dataset.")
    eigenfuncs = np.array([{}] * len(ef_idxs), dtype=dict)
    with open(dataset.datfile, "rb") as istream:
        for dict_idx, ef_idx in enumerate(ef_idxs):
            efs = read_eigenfunction(istream, dataset.header, ef_idx)
            efs.update({"eigenvalue": dataset.eigenvalues[ef_idx]})
            eigenfuncs[dict_idx] = efs
    return eigenfuncs


class EigenfunctionHandler:
    def __init__(self, data, ef_ax):
        self.data = data
        self.ef_ax = ef_ax
        # check presence of eigenfunctions
        error_msg = "None of the given datfiles has eigenfunctions written to it."
        if isinstance(self.data, LegolasDataSeries):
            if not any(self.data.efs_written):
                raise EigenfunctionsNotPresent(error_msg)
            self._mark_eigenfunctions_present()
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

    def on_point_pick(self, event):
        artist = event.artist
        # if artist is a legend item, return (this attribute has been set manually)
        if hasattr(artist, "is_legend_item"):
            return
        # retrieve figure and axis for this artist, save limits to prevent zoom reset
        fig, ax = artist.figure, artist.axes
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        # Retrieves the dataset associated with the current points. This has been
        # set when the spectrum was added to the plot.
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
            self._selected_idxs.update({associated_ds: {f"{idx}": marked_point}})
            # checks if eigenfunctions are available
            self._check_eigenfunction_availability()
        # handle right clicking
        elif event.mouseevent.button == 3:
            # remove selected index from list
            selected_artist = self._selected_idxs.get(associated_ds, {}).pop(
                str(idx), None
            )
            if selected_artist is not None:
                selected_artist.remove()
        # this fixes a zoom reset when picking the figure
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        fig.canvas.draw()

    def on_key_press(self, event):
        # if nothing is selected, return
        if not self._selected_idxs:
            return
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
        self.ef_ax.clear()
        ef_name = self._ef_names[self._selected_ef_name_idx]
        for ds, idxs_dict in self._selected_idxs.items():
            idxs = np.array([int(idx) for idx in idxs_dict.keys()])
            ef_container = get_eigenfunctions(ds, idxs)

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
        self._selected_ef_name_idx += 1
        if self._selected_ef_name_idx > len(self._ef_names) - 1:
            self._selected_ef_name_idx = 0

    def _select_previous_eigenfunction(self):
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
        becomes :math:`(x, y, z)` or :math:`(r, \theta, z)`.

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

    def _mark_eigenfunctions_present(self):
        return
        raise NotImplementedError

    def _check_eigenfunction_availability(self):
        """
        Checks if eigenfunctions are available. This is important for
        multispectra, where every dataset does not necessarily has eigenfunctions
        associated with it. This routine will check the marked points, if a point
        is encountered without eigenfunctions it will be removed and a warning
        will be thrown.
        """
        return
        if isinstance(self.data, LegolasDataSeries):
            raise NotImplementedError
