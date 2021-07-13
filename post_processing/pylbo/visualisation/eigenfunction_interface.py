import abc
import numpy as np
from matplotlib.collections import PathCollection
from pylbo.data_containers import LegolasDataSeries
from pylbo.utilities.toolbox import add_pickradius_to_item


class EigenfunctionInterface:
    __metaclass__ = abc.ABCMeta

    def __init__(self, data, axis):
        self.data = data
        self.axis = axis
        self._check_data_is_present()
        # holds the points that are currently selected in the form of a "double" dict:
        # {"ds instance" : {"index" : line2D instance}}
        self._selected_idxs = {}
        self._use_real_part = True
        self._selected_name_idx = 0
        self._function_names = None
        self._retransform = False

        self._condition_to_make_transparent = None
        self._transparent_data = False
        self._unmarked_alpha = None

    def _check_data_is_present(self):
        """
        Checks if the required data is present to draw for example
        eigenfunctions, is overloaded in subclasses.
        """
        pass

    def _artist_has_valid_attributes(self, event):
        """
        Checks if a given event has valid attributes, this prevents triggering the
        interface when clicking on legend items, for example.

        Parameters
        ----------
        event : ~matplotlib.backend_bases.PickEvent
            The pick event.

        Returns
        -------
        bool
            `True` if all conditions are met and callbacks can be connected, `False`
            otherwise.
        """
        # this attr has been set in case of multiple axes, prevents double event trigger
        if hasattr(self, "associated_ds_ax"):
            if not event.mouseevent.inaxes == self.associated_ds_ax:
                return False
        artist = event.artist
        # this attr has been set when adding the legend, return for legend items
        if hasattr(artist, "is_legend_item"):
            return False
        if not hasattr(artist, "dataset"):
            return False
        return True

    def _clear_figure_and_selection(self):
        """
        Clears the current figure, clears the dictionary of selected eigenvalues.
        """
        for ds_artist_dict in self._selected_idxs.values():
            for artist in ds_artist_dict.values():
                artist.remove()
        self._selected_idxs.clear()
        self.axis.clear()

    def _switch_real_and_imaginary_part(self):
        """
        Switches between the real and imaginary part of a given function.
        """
        self._use_real_part = not self._use_real_part

    def _select_next_function(self):
        """
        Increments the index of the currently selected function by 1.
        """
        self._selected_name_idx += 1
        if self._selected_name_idx > len(self._function_names) - 1:
            self._selected_name_idx = 0

    def _select_previous_function(self):
        """
        Decrements the index of the currently selected function by 1.
        """
        self._selected_name_idx -= 1
        if self._selected_name_idx < 0:
            self._selected_name_idx = len(self._function_names) - 1

    def _retransform_functions(self):
        """
        Toggles a retransform of a function, for example an eigenfunction
        :math:`v_r \\leftrightarrow rv_r`"
        """
        self._retransform = not self._retransform

    def _print_selected_eigenvalues(self):
        """
        Prints all selected eigenvalues to the console as an array.
        """
        if not self._selected_idxs:
            return
        print("Currently selected eigenvalues:")
        for ds, points in self._selected_idxs.items():
            idxs = np.array([int(idx) for idx in points.keys()])
            print(f"{ds.datfile.stem} | {ds.eigenvalues[idxs]}")

    def _get_label(self, ds, ev_idx, w):
        """
        Returns the label used in the legend. In case of a data series, the datfile
        name is prepended.

        Parameters
        ----------
        ds : ~pylbo.data_containers.LegolasDataSet
            The current dataset
        ev_idx : int
            The index of the current eigenvalue in the corresponding array
        w : float, complex
            The eigenvalue to use in the label

        Returns
        -------
        str
            The label to use in the legend.
        """
        label = rf"$\omega_{{{ev_idx}}}$ = {w:2.3e}"
        if isinstance(self.data, LegolasDataSeries):
            label = f"{ds.datfile.stem} | {label}"
        return label

    @abc.abstractmethod
    def _get_title(self):
        """
        Creates the title of a given plot, has to be overridden in a subclass.
        """
        pass

    @abc.abstractmethod
    def update_plot(self):
        """
        Updates the plot when an event is triggered, clears and then redraws
        the functions. Rescaling of the axes is done automatically.
        Has to be overridden in a subclass.
        """
        pass

    def on_point_pick(self, event):
        """
        Determines what happens when an eigenvalue is clicked.

        Parameters
        ----------
        event : ~matplotlib.backend_bases.PickEvent
            The pick event.
        """
        if not self._artist_has_valid_attributes(event):
            return
        # retrieve limits to prevent resetting zoom
        ax = event.artist.axes
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        if event.mouseevent.button == 1:
            self.on_left_click(event)
        elif event.mouseevent.button == 3:
            self.on_right_click(event)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        event.artist.figure.canvas.draw()

    def on_key_press(self, event):
        """
        Determines what happens when a key is pressed.

        Parameters
        ----------
        event : ~matplotlib.backend_bases.KeyEvent
            The key event.
        """
        if event.key == "m":
            self._mark_points_without_data_written()
        elif event.key == "d":
            self._clear_figure_and_selection()
        elif event.key == "i":
            self._switch_real_and_imaginary_part()
        elif event.key == "up":
            self._select_next_function()
        elif event.key == "down":
            self._select_previous_function()
        elif event.key == "t":
            self._retransform_functions()
        elif event.key == "w":
            self._print_selected_eigenvalues()
        if event.key in ("enter", "up", "down", "i", "t"):
            self.update_plot()
        event.canvas.draw()

    def on_left_click(self, event):
        """
        Determines what happens when left-clicking an eigenvalue.

        Parameters
        ----------
        event : ~matplotlib.backend_bases.PickEvent
            The pick event.
        """
        idx, xdata, ydata = self._get_clicked_point_data(event)
        associated_ds = event.artist.dataset
        # skip if point index is already in list
        if str(idx) in self._selected_idxs.get(associated_ds, {}).keys():
            return
        (marked_point,) = event.artist.axes.plot(
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

    def on_right_click(self, event):
        """
        Determines what happens when right-clicking an eigenvalue.

        Parameters
        ----------
        event : ~matplotlib.backend_bases.PickEvent
            The pick event.
        """
        idx, _, _ = self._get_clicked_point_data(event)
        # remove selected index from list
        associated_ds = event.artist.dataset
        selected_artist = self._selected_idxs.get(associated_ds, {}).pop(str(idx), None)
        if selected_artist is not None:
            selected_artist.remove()
            # if no items remaining for this ds, remove key
            if len(self._selected_idxs[associated_ds]) == 0:
                self._selected_idxs.pop(associated_ds)

    def _get_clicked_point_data(self, event):
        """
        Retrieves the index (in the eigenvalue array), x data coordinate and y data
        coordinate of the eigenvalue nearest to the clicked point.

        Parameters
        ----------
        event : ~matplotlib.backend_bases.PickEvent
            The pick event.

        Returns
        -------
        idx : integer
            The index of the selected point in the eigenvalue array
        xdata : float
            The x data coordinate of the selected eigenvalue
        ydata : float
            The y data coordinate of the selected eigenvalue
        """
        artist = event.artist
        idxs = event.ind
        if isinstance(artist, PathCollection):
            # this block is entered for points drawn using scatter instead of plot
            xdata, ydata = np.split(artist.get_offsets(), [-1], axis=1)
            xdata = xdata.squeeze()
            ydata = ydata.squeeze()
        else:
            xdata = artist.get_xdata()
            ydata = artist.get_ydata()
        # In case of overlap between clicked points then multiple
        # indices were selected, so we check smallest distance to mouse click
        if len(idxs) == 1:
            idx = idxs[0]
        else:
            mouse_x = event.mouseevent.xdata
            mouse_y = event.mouseevent.ydata
            distances = (mouse_x - xdata[idxs]) ** 2 + (mouse_y - ydata[idxs]) ** 2
            idx = idxs[distances.argmin()]
        return idx, xdata[idx], ydata[idx]

    def _mark_points_without_data_written(self):
        """
        For dataseries, it is possible that not all datasets in the series
        have eigenfunctions associated with them. This routine will toggle a change
        in the opacity value for datapoints with no functions, so they are
        clearly distinguishable from those who do have them.
        """
        if not isinstance(self.data, LegolasDataSeries):
            return
        if all(getattr(self.data, self._condition_to_make_transparent)):
            return
        self._transparent_data = not self._transparent_data
        for ax in self.axis.figure.get_axes():
            for child in ax.get_children():
                # the ones with this attribute are all vertical rows of datapoints
                if hasattr(child, "dataset"):
                    if self._unmarked_alpha is None:
                        self._unmarked_alpha = child.get_alpha()
                    if not getattr(child.dataset, self._condition_to_make_transparent):
                        child.set_alpha(0)
                    else:
                        child.set_alpha(self._unmarked_alpha)
