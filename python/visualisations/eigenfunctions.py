from utilities.plot_data import plot_eigenfunctions

import numpy as np
import matplotlib.pyplot as plt

class EigenfunctionHandler:
    def __init__(self, data, fig, ax):
        print("")
        print("-" * 50)
        print(">>> INTERACTIVE PLOTTING OF EIGENFUNCTIONS <<<")
        print("- Left-click  : select spectrum points (red cross)")
        print("- Right-click : deselect spectrum points")
        print("- Enter       : plot eigenfunctions corresponding to selected points")
        print("- Up-arrow    : cycle upwards through eigenfunction variables")
        print("- Down-arrow  : cycle downwards through eigenfunction variables")
        print("- d-key       : clears selection and eigenfunction figure")
        print("- x-key       : close figures and finish program ")
        print("-" * 50)
        print("")

        self.spec_fig = fig
        self.spec_ax = ax
        self.omegas = data.omegas

        self.w_idx_list = []
        self.current_var = 0
        self.plot_real = True

        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))

        self.ef_list = data.ef_list
        self.grid = data.grid
        self.eigenfuncs = data.eigenfuncs


    def connect_interactive_funcs(self):
        def _on_clicking(event):
            toolbar = self.spec_fig.canvas.manager.toolbar
            # Toolbar should be deactivated when in interactive mode
            if toolbar.mode == '':
                if event.button == 1:
                    # Prevent error when clicking outside of axis window
                    if event.xdata is None or event.ydata is None:
                        return
                    # Search nearest spectrum point to click, then plot
                    idx = self._find_spectrum_point(event.xdata, event.ydata)
                    # Check for double plotting of points
                    if idx not in self.w_idx_list:
                        self.w_idx_list.append(idx)
                        self.spec_ax.plot(np.real(self.omegas[idx]), np.imag(self.omegas[idx]), 'rx',
                                          markersize=8, picker=10, label='w_point')
            self.spec_fig.canvas.draw()

        def _on_picking(event):
            # Right clicking on point removes the point. Radius in which clicking is detected depends on 'picker'
            # argument when plotting the point itself (see _on_clicking)
            if event.mouseevent.button == 3:
                if hasattr(event.artist, 'get_label') and event.artist.get_label() == 'w_point':
                    event.artist.remove()
                    # remove point from index list as well
                    idx = self._find_spectrum_point(event.artist.get_xdata(), event.artist.get_ydata())
                    self.w_idx_list.remove(idx)
            self.spec_fig.canvas.draw()

        def _on_typing(event):
            # Press 'x' to close figures
            if event.key == 'x':
                for i in plt.get_fignums():
                    plt.close(i)
                return

            # Press 'd' clears current figure and selection
            if event.key == 'd':
                self.w_idx_list.clear()
                for artist in self.spec_ax.get_children():
                    if hasattr(artist, 'get_label') and artist.get_label() == 'w_point':
                        artist.remove()
                self.ax.clear()
                self.fig.set_visible(False)
                self.fig.canvas.draw()
                self.spec_fig.canvas.draw()

            # Do nothing if nothing is selected
            if not self.w_idx_list:
                return

            # Press 'i' to switch between real and imaginary parts
            if event.key == 'i':
                self._switch_real()
                self.ax.clear()
                var = self._get_current_variable()
                plot_eigenfunctions(self.fig, self.ax, self.omegas, self.grid, self.eigenfuncs,
                                    var, self.w_idx_list, self.plot_real)
            if event.key == 'down':
                self.var = self._select_prev_variable()
            if event.key == 'up':
                self.var = self._select_next_variable()

            # Press 'enter' to plot eigenfunctions, press up/down to cycle
            if event.key == 'enter' or event.key == 'up' or event.key == 'down' or event.key == 'i':
                self.fig.set_visible(True)
                self.ax.clear()
                var = self._get_current_variable()
                plot_eigenfunctions(self.fig, self.ax, self.omegas, self.grid, self.eigenfuncs,
                                    var, self.w_idx_list, self.plot_real)
            self.fig.canvas.draw()
            self.spec_fig.canvas.draw()

        # connect interactive functions
        self.fig.canvas.mpl_connect('key_press_event', _on_typing)
        self.spec_fig.canvas.mpl_connect('key_press_event', _on_typing)
        self.spec_fig.canvas.mpl_connect('button_press_event', _on_clicking)
        self.spec_fig.canvas.mpl_connect('pick_event', _on_picking)
        self.fig.set_visible(False)


    def _select_next_variable(self):
        self.current_var += 1
        if self.current_var > 7:
            self.current_var = 0
        return self.ef_list[self.current_var]


    def _select_prev_variable(self):
        self.current_var -= 1
        if self.current_var < 0:
            self.current_var = 7
        return self.ef_list[self.current_var]


    def _get_current_variable(self):
        return self.ef_list[self.current_var]


    def _switch_real(self):
        self.plot_real = not self.plot_real


    def _find_spectrum_point(self, x, y):
        w_clicked = x + y*1j
        idx = (np.abs(self.omegas - w_clicked)).argmin()
        return idx



