import matplotlib.pyplot as plt
import matplotlib.lines as mpl_lines
import matplotlib.patches as mpl_patches
import numpy as np

from data_container import LEGOLASDataContainer
from visualisations.eigenfunctions import EigenfunctionHandler
from utilities.plot_data import plot_spectrum, plot_matrix
from utilities.continua import get_continuum_regions


class SingleSpectrum(LEGOLASDataContainer):
    def __init__(self, namelist_array):
        super().__init__(namelist_array)

        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))
        self.data = self.datacontainer.pop()

    def plot(self):
        spectrum_title = self.data.current_eq + '   ({}, {} gridpts)'.format(self.data.geometry, self.data.gridpts)
        plot_spectrum(self.fig, self.ax, self.data.omegas, title=spectrum_title)

        if self.data.equil_data is not None:
            self._annotate_frequencies()
        if self.data.show_mats:
            self._plot_matrices()
        if self.data.show_equil:
            self._plot_equilibria()
        if self.data.show_eigenfuncs:
            efh = EigenfunctionHandler(self.data, self.fig, self.ax)
            efh.connect_interactive_funcs()

    def _plot_matrices(self):
        fig, ax = plt.subplots(1, 2, figsize=(12, 8))
        plot_matrix(fig, ax, 0, self.data.matA, title='Matrix A')
        plot_matrix(fig, ax, 1, self.data.matB, title='Matrix B')

    def _plot_equilibria(self):
        fig, ax = plt.subplots(1, figsize=(12, 8))
        colors = ['blue', 'red', 'green', 'cyan', 'orange', 'grey']

        for idx, var in enumerate(['rho0', 'T0', 'v02', 'v03', 'B02', 'B03']):
            ax.plot(self.data.grid_gauss, self.data.equil_data[var], 'o', markersize=2, color=colors[idx], label=var)
            ax.plot(self.data.grid_gauss, self.data.equil_data[var], color=colors[idx], alpha=0.3)
        ax.set_xlabel('Gaussian grid')
        ax.set_title('Equilibrium configuration')
        ax.set_xlim([self.data.grid_gauss[0], self.data.grid_gauss[-1]])
        ax.legend(loc='best')

    def _annotate_frequencies(self):
        legend_items = []
        alpha_box = 0.3
        alpha_point = 0.5

        def on_legend_pick(event):
            artist = event.artist
            # do nothing if clicking outside of legend
            if not artist in regions:
                return
            item = regions[artist]
            visible = not item.get_visible()
            item.set_visible(visible)
            if visible:
                if isinstance(artist, mpl_lines.Line2D):
                    artist.set_alpha(alpha_point)
                else:
                    artist.set_alpha(alpha_box)
            else:
                artist.set_alpha(0.1)
            self.fig.canvas.draw()

        def draw_region(array, facecolor, legend_lbl=None, draw_horizontal=False):
            lb = np.min(array)
            rb = np.max(array)
            if np.abs(lb - rb) > 1e-10:
                if draw_horizontal:
                    legend_item = self.ax.axhspan(np.min(array), np.max(array), facecolor=facecolor,
                                                  alpha=alpha_box, label=legend_lbl)
                else:
                    legend_item = self.ax.axvspan(np.min(array), np.max(array), facecolor=facecolor,
                                                  alpha=alpha_box, label=legend_lbl)
            else:   # handles the case where continua are collapsed to a single point
                legend_item, = self.ax.plot([lb], 0, marker='p', markersize=8, color=facecolor,
                                            alpha=alpha_point, label=legend_lbl)
            legend_items.append(legend_item)

        # calculate all continuum regions
        wS_pos, wS_neg, wA_pos, wA_neg, wTH = get_continuum_regions(self.data)
        # draw regions on plot
        draw_region(wS_pos, facecolor='red', legend_lbl=r'$\Omega_S^+$ slow continuum')
        draw_region(wS_neg, facecolor='red', legend_lbl=r'$\Omega_S^-$ slow continuum')
        draw_region(wA_pos, facecolor='blue', legend_lbl=r'$\Omega_A^+$ Alfvén continuum')
        draw_region(wA_neg, facecolor='blue', legend_lbl=r'$\Omega_A^-$ Alfvén continuum')
        # thermal continuum is imaginary, so draw horizontally
        draw_region(wTH, facecolor='green', legend_lbl=r'$\Omega_{th}$ thermal continuum', draw_horizontal=True)
        legend = self.ax.legend(loc='best')

        # retrieve legend patches and items for interactive enabling of continuum regions
        regions = {}
        patches = legend.get_patches()
        lines = legend.get_lines()
        for region in legend_items:
            if isinstance(region, mpl_patches.Polygon):
                l_item = patches.pop(0)
            elif isinstance(region, mpl_lines.Line2D):
                l_item = lines.pop(0)
            else:
                raise ValueError
            l_item.set_picker(5)
            region.set_visible(False)   # make continuum regions invisible by default
            regions[l_item] = region
        self.fig.canvas.mpl_connect('pick_event', on_legend_pick)
        self._plot_continua(wS_pos, wS_neg, wA_pos, wA_neg, wTH)

    def _plot_continua(self, wS_pos, wS_neg, wA_pos, wA_neg, wTH):
        # plot every continuum on a separate figure
        fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex='all')

        ax[0].plot(self.data.grid_gauss, wS_pos, color='red', linestyle='solid',
                   label=r'$\Omega_S^-$ slow continuum')
        ax[0].plot(self.data.grid_gauss, wS_neg, color='red', linestyle='dashed',
                   label=r'$\Omega_S^+$ slow continuum')
        ax[1].plot(self.data.grid_gauss, wA_pos, color='blue', linestyle='solid',
                   label=r'$\Omega_A^-$ Alfvén continuum')
        ax[1].plot(self.data.grid_gauss, wA_neg, color='blue', linestyle='dashed',
                   label=r'$\Omega_A^+$ Alfvén continuum')
        ax[2].plot(self.data.grid_gauss, wTH, color='green', linestyle='solid',
                   label=r'Im$(\Omega_{th})$ thermal continuum')
        ax[2].set_xlabel('Gaussian grid')
        for ax_panel in ax:
            ax_panel.legend(loc='best')
            ax_panel.set_xlim([self.data.grid_gauss[0], self.data.grid_gauss[-1]])
        fig.subplots_adjust(hspace=0.05)
