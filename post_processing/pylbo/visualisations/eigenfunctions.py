import matplotlib.pyplot as plt
import numpy as np

class EigenfunctionHandler:
    def __init__(self, ps, merge_figs):
        self.ps = ps
        self.ef_grid, self.eigenfunctions = self.ps.ds.get_eigenfunctions()
        if merge_figs:
            self.fig, ax = plt.subplots(1, 2, figsize=(14, 6))
            # remove 'old' spectrum figure
            plt.delaxes(self.ps.ax)
            plt.close(self.ps.fig)
            self.ax = ax[0]
            self.ps_ax = ax[1]
            self.ps_fig = self.fig
            # replot spectrum, now using the subplot
            self.ps._plot(self.fig, self.ps_ax)
        else:
            self.fig, self.ax = plt.subplots(1, figsize=(12, 6))
            self.ps_ax = self.ps.ax
            self.ps_fig = self.ps.fig
        self.ef_names = self.ps.ds.header.get('ef_names')
        self.ev_indices = []
        self.selected_ef_idx = 0
        self.ef_real_part = True

    def connect_figure_events(self):
        toolbar = self.ps_fig.canvas.manager.toolbar
        def on_left_click(event):
            # do nothing if toolbar is activated
            if not toolbar.mode == '':
                return
            if event.button == 1:
                # do nothing if clicking outside of window
                if event.xdata is None or event.ydata is None:
                    return
                # search nearest spectrum point to click
                index = self._find_spectrum_point_index(event.xdata, event.ydata)
                if index is None:
                    return
                # only plot if point is not already selected
                if index not in self.ev_indices:
                    self.ev_indices.append(index)
                    self.ps_ax.plot(self.ps.ds.eigenvals[index].real, self.ps.ds.eigenvals[index].imag,
                                    'rx', markersize=8, picker=10, label='eigenvalue')
            self.ps_fig.canvas.draw()

        def on_right_click(event):
            if not toolbar.mode == '':
                return
            if event.mouseevent.button == 3:
                if hasattr(event.artist, 'get_label') and event.artist.get_label() == 'eigenvalue':
                    # right clicking returns (x, y) as ([x], [y]) so use unpacking
                    index = self._find_spectrum_point_index(*event.artist.get_xdata(), *event.artist.get_ydata())
                    if index is None:
                        return
                    self.ev_indices.remove(index)
                    event.artist.remove()
            self.ps_fig.canvas.draw()

        def on_key_press(event):
            # 'x' closes figure
            if event.key == 'x':
                for i in plt.get_fignums():
                    plt.close(i)
                return
            # 'd' clears figure and selection
            if event.key == 'd':
                self.ev_indices.clear()
                for artist in self.ps_ax.get_children():
                    if hasattr(artist, 'get_label') and artist.get_label() == 'eigenvalue':
                        artist.remove()
                self.ax.clear()
                self.fig.canvas.draw()
                self.ps_fig.canvas.draw()
            # if nothing is selected, return
            if not self.ev_indices:
                return
            # 'i' switches between real and imaginary parts
            if event.key == 'i':
                self.ef_real_part = not self.ef_real_part
            # 'up'/'down' cycles between eigenfunctions
            if event.key == 'up':
                self._select_next_ef_index()
            if event.key == 'down':
                self._select_prev_ef_index()
            if event.key in ('enter', 'up', 'down', 'i'):
                self._plot_eigenfunctions()

        self.ps_fig.canvas.mpl_connect('button_press_event', on_left_click)
        self.ps_fig.canvas.mpl_connect('pick_event', on_right_click)
        self.ps_fig.canvas.mpl_connect('key_press_event', on_key_press)
        self.fig.canvas.mpl_connect('key_press_event', on_key_press)

    def _find_spectrum_point_index(self, x, y):
        # get distance from (x, y) to all points
        distances = np.sqrt((self.ps.ds.eigenvals.real - x)**2 + (self.ps.ds.eigenvals.imag - y)**2)
        # index of point with closest distance
        idx = distances.argmin()
        ev_found = self.ps.ds.eigenvals[idx]
        # calculate (x, y) of point in pixels. (0, 0) is bottom-left of figure
        ev_x_pixels, ev_y_pixels = self.ps_ax.transData.transform((ev_found.real, ev_found.imag))
        click_x_pixels, click_y_pixels = self.ps_ax.transData.transform((x, y))
        # only select point if clicked within certain distance, say 15 pixels
        pixel_criterion = 15
        dist_pixels = np.sqrt((ev_x_pixels - click_x_pixels)**2 + (ev_y_pixels - click_y_pixels)**2)
        if dist_pixels < pixel_criterion:
            return idx
        return None

    def _plot_eigenfunctions(self):
        self.ax.clear()
        ef_name = self.ef_names[self.selected_ef_idx]
        for ev_idx in self.ev_indices:
            # eigenfunctions are contained in columns
            ef = self.eigenfunctions.get(ef_name)[:, ev_idx]
            if self.ef_real_part:
                ef = ef.real
            else:
                ef = ef.imag
            label = r'$\omega${} = {:2.5e}'.format(ev_idx, self.ps.ds.eigenvals[ev_idx])
            self.ax.plot(self.ef_grid, ef, label=label)
        self.ax.axhline(y=0, linestyle='dotted', color='grey')
        self.ax.axvline(x=0, linestyle='dotted', color='grey')
        if self.ef_real_part:
            title = '{} eigenfunction (real part)'.format(ef_name)
        else:
            title = '{} eigenfunction (imag part)'.format(ef_name)
        self.ax.set_title(title)
        self.ax.legend(loc='best')
        self.fig.canvas.draw()

    def _select_next_ef_index(self):
        self.selected_ef_idx += 1
        if self.selected_ef_idx > len(self.ef_names) - 1:
            self.selected_ef_idx = 0

    def _select_prev_ef_index(self):
        self.selected_ef_idx -= 1
        if self.selected_ef_idx < 0:
            self.selected_ef_idx = len(self.ef_names) - 1
