from data_container import LEGOLASDataContainer
from visualisations.eigenfunctions import EigenfunctionHandler
from utilities.plot_data import plot_spectrum, plot_matrix

import matplotlib.pyplot as plt
import matplotlib.lines as mpl_lines
import matplotlib.patches as mpl_patches
import numpy as np


class SingleSpectrum(LEGOLASDataContainer):
    def __init__(self, namelist_array):
        super().__init__(namelist_array)

        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))
        self.data = self.datacontainer.pop()

    def plot(self):
        spectrum_title = self.data.current_eq + '   ({})'.format(self.data.geometry)
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

        for var in ['rho0', 'T0', 'v02', 'v03', 'B02', 'B03']:
            ax.plot(self.data.grid_gauss, self.data.equil_data[var], label=var)
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

        def draw_region(array, facecolor, legend_lbl=None):
            lb = np.min(array)
            rb = np.max(array)

            if np.abs(lb - rb) > 1e-10:
                legend_item = self.ax.axvspan(np.min(array), np.max(array), facecolor=facecolor,
                                              alpha=alpha_box, label=legend_lbl)
            else:   # handles the case where continua are collapsed to a single point
                legend_item, = self.ax.plot([lb], 0, marker='p', markersize=8, color=facecolor,
                                            alpha=alpha_point, label=legend_lbl)

            legend_items.append(legend_item)

        k2 = self.data.params['k2']
        k3 = self.data.params['k3']

        rho = self.data.equil_data['rho0']
        B02 = self.data.equil_data['B02']
        B03 = self.data.equil_data['B03']
        B0 = np.sqrt(B02**2 + B03**2)
        v02 = self.data.equil_data['v02']
        v03 = self.data.equil_data['v03']

        gamma = self.data.gamma
        p = rho * self.data.equil_data['T0']

        # Omega_0 = k0 dot v (doppler shift)
        doppler_shift = k2 * v02 + k3 * v03

        if self.data.geometry == 'Cartesian':
            wA2 = (1/rho) * (k2*B02 + k3*B03)**2
            wS2 = (gamma*p / (gamma*p + B0**2)) * wA2
        else:
            r = self.data.grid_gauss
            F = (k2 * B02 / r) + k3 * B03
            wA2 = F**2 / rho
            wS2 = (gamma*p / (gamma*p + B0**2)) * F**2 / rho

        alfven_cont_neg = doppler_shift - np.sqrt(wA2)
        alfven_cont_pos = doppler_shift + np.sqrt(wA2)
        slow_cont_neg = doppler_shift - np.sqrt(wS2)
        slow_cont_pos = doppler_shift + np.sqrt(wS2)

        draw_region(alfven_cont_neg, facecolor='red', legend_lbl=r'$\Omega_A^-$ Alfvén continuum')
        draw_region(alfven_cont_pos, facecolor='green', legend_lbl=r'$\Omega_A^+$ Alfvén continuum')
        draw_region(slow_cont_neg, facecolor='blue', legend_lbl=r'$\Omega_S^-$ slow continuum')
        draw_region(slow_cont_pos, facecolor='yellow', legend_lbl=r'$\Omega_S^+$ slow continuum')
        legend = self.ax.legend(loc='best')

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
        self._plot_continua(alfven_cont_neg, alfven_cont_pos, slow_cont_neg, slow_cont_pos)

    def _plot_continua(self, alfven_cont_neg, alfven_cont_pos, slow_cont_neg, slow_cont_pos):
        # plot Alfven and slow continua on a separate figure
        fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex='all')

        ax[0].plot(self.data.grid_gauss, alfven_cont_neg, color='red', label=r'$\Omega_A^-$ Alfvén continuum')
        ax[0].plot(self.data.grid_gauss, alfven_cont_pos, color='green', label=r'$\Omega_A^+$ Alfvén continuum')
        ax[0].legend(loc='best')
        ax[0].set_xlim([self.data.grid_gauss[0], self.data.grid_gauss[-1]])

        ax[1].plot(self.data.grid_gauss, slow_cont_neg, color='blue', label=r'$\Omega_S^-$ slow continuum')
        ax[1].plot(self.data.grid_gauss, slow_cont_pos, color='yellow', label=r'$\Omega_S^+$ slow continuum')
        ax[1].set_xlabel('Gaussian grid')
        ax[1].legend(loc='best')
        ax[1].set_xlim([self.data.grid_gauss[0], self.data.grid_gauss[-1]])
        fig.subplots_adjust(hspace=0)
