import matplotlib.pyplot as plt
import numpy as np

class EigenfunctionHandler:
    def __init__(self, ds):
        self.ds = ds
        self.ef_grid, self.eigenfunctions = self.ds.get_eigenfunctions()
        self.ef_names = self.ds.header.get('ef_names')
        self._ev_indices = []
        self._selected_ef_idx = 0
        self._ef_real_part = True
        self._ef_fig = None
        self._ef_ax = None

    def connect_figure_events(self, spectrum, ef_fig, ef_ax):
        toolbar = spectrum.fig.canvas.manager.toolbar
        self._ef_fig = ef_fig
        self._ef_ax = ef_ax

        def on_left_click(event):
            # do nothing if toolbar is activated
            if not toolbar.mode == '':
                return
            if event.button == 1:
                # do nothing if clicking outside of window
                if event.xdata is None or event.ydata is None:
                    return
                # search nearest spectrum point to click
                index = self._find_spectrum_point_index(spectrum, event.xdata, event.ydata)
                if index is None:
                    return
                # only plot if point is not already selected
                if index not in self._ev_indices:
                    self._ev_indices.append(index)
                    spectrum.ax.plot(self.ds.eigenvalues[index].real, self.ds.eigenvalues[index].imag,
                                     'rx', markersize=8, picker=10, label='eigenvalue')
            spectrum.fig.canvas.draw()

        def on_right_click(event):
            if not toolbar.mode == '':
                return
            if event.mouseevent.button == 3:
                if hasattr(event.artist, 'get_label') and event.artist.get_label() == 'eigenvalue':
                    # right clicking returns (x, y) as ([x], [y]) so use unpacking
                    index = self._find_spectrum_point_index(spectrum, *event.artist.get_xdata(),
                                                            *event.artist.get_ydata())
                    if index is None:
                        return
                    self._ev_indices.remove(index)
                    event.artist.remove()
            spectrum.fig.canvas.draw()

        def on_key_press(event):
            # 'x' closes figure
            if event.key == 'x':
                for i in plt.get_fignums():
                    plt.close(i)
                return
            # 'd' clears figure and selection
            if event.key == 'd':
                self._ev_indices.clear()
                for artist in spectrum.ax.get_children():
                    if hasattr(artist, 'get_label') and artist.get_label() == 'eigenvalue':
                        artist.remove()
                self._ef_ax.clear()
                self._ef_fig.canvas.draw()
                spectrum.fig.canvas.draw()
            # if nothing is selected, return
            if not self._ev_indices:
                return
            # 'i' switches between real and imaginary parts
            if event.key == 'i':
                self._ef_real_part = not self._ef_real_part
            # 'up'/'down' cycles between eigenfunctions
            if event.key == 'up':
                self._select_next_ef_index()
            if event.key == 'down':
                self._select_prev_ef_index()
            if event.key in ('enter', 'up', 'down', 'i'):
                self._plot_eigenfunctions()

        spectrum.fig.canvas.mpl_connect('button_press_event', on_left_click)
        spectrum.fig.canvas.mpl_connect('pick_event', on_right_click)
        spectrum.fig.canvas.mpl_connect('key_press_event', on_key_press)
        self._ef_fig.canvas.mpl_connect('key_press_event', on_key_press)

    def _find_spectrum_point_index(self, spectrum, x, y):
        # this will return a single index/eigenvalue, since we supply 1 'guess' (x, y)
        idx, ev_found = self.ds.get_nearest_eigenvalues(complex(x, y))
        idx, ev_found = *idx, *ev_found
        # calculate (x, y) of point in pixels. (0, 0) is bottom-left of figure
        ev_x_pixels, ev_y_pixels = spectrum.ax.transData.transform((ev_found.real, ev_found.imag))
        click_x_pixels, click_y_pixels = spectrum.ax.transData.transform((x, y))
        # only select point if clicked within certain distance, say 15 pixels
        pixel_criterion = 15
        dist_pixels = np.sqrt((ev_x_pixels - click_x_pixels)**2 + (ev_y_pixels - click_y_pixels)**2)
        if dist_pixels < pixel_criterion:
            return idx
        return None

    def _plot_eigenfunctions(self):
        self._ef_ax.clear()
        ef_name = self.ef_names[self._selected_ef_idx]
        for ev_idx in self._ev_indices:
            # eigenfunctions are contained in columns
            ef = self.eigenfunctions.get(ef_name)[:, ev_idx]
            if self._ef_real_part:
                ef = ef.real
            else:
                ef = ef.imag
            label = r'$\omega${} = {:2.5e}'.format(ev_idx, self.ds.eigenvalues[ev_idx])
            self._ef_ax.plot(self.ef_grid, ef, label=label)
        self._ef_ax.axhline(y=0, linestyle='dotted', color='grey')
        self._ef_ax.axvline(x=self.ds.x_start, linestyle='dotted', color='grey')
        if self._ef_real_part:
            title = '{} eigenfunction (real part)'.format(ef_name)
        else:
            title = '{} eigenfunction (imag part)'.format(ef_name)
        self._ef_ax.set_title(title)
        self._ef_ax.legend(loc='best')
        self._ef_fig.canvas.draw()

    def _select_next_ef_index(self):
        self._selected_ef_idx += 1
        if self._selected_ef_idx > len(self.ef_names) - 1:
            self._selected_ef_idx = 0

    def _select_prev_ef_index(self):
        self._selected_ef_idx -= 1
        if self._selected_ef_idx < 0:
            self._selected_ef_idx = len(self.ef_names) - 1

    def retrieve_eigenfunctions(self, ev_guesses):
        result = {ef_name: [] for ef_name in self.ef_names}
        result.update({'eigenvals': []})
        idxs, evs = self.ds.get_nearest_eigenvalues(ev_guesses)
        for idx, ev in zip(idxs, evs):
            result['eigenvals'].append(ev)
            for ef_name in self.ef_names:
                result.get(ef_name).append(self.eigenfunctions.get(ef_name)[:, idx])
        return result






