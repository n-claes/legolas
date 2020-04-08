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
        self.eq_fig, self.eq_ax = plt.subplots(1, figsize=(12, 8))
        self.co_fig, self.co_ax = plt.subplots(3, 1, figsize=(8, 8), sharex='all')

        self.data = self.datacontainer.pop()

    def plot(self):
        spectrum_title = self.data.current_eq + '   ({}, {} gridpts)'.format(self.data.geometry, self.data.gridpts)
        plot_spectrum(self.fig, self.ax, self.data.omegas, title=spectrum_title)

        if self.data.equil_data is not None:
            self._annotate_frequencies()
        else:
            self.co_fig.close()
        if self.data.show_mats:
            self._plot_matrices()
        if self.data.show_equil:
            self._plot_equilibria()
        else:
            self.eq_fig.close()
        if self.data.show_eigenfuncs:
            efh = EigenfunctionHandler(self.data, self.fig, self.ax)
            efh.connect_interactive_funcs()

    def _plot_matrices(self):
        fig, ax = plt.subplots(1, 2, figsize=(12, 8))
        plot_matrix(fig, ax, 0, self.data.matA, title='Matrix A')
        plot_matrix(fig, ax, 1, self.data.matB, title='Matrix B')

    def _plot_equilibria(self):
        colors = ['blue', 'red', 'green', 'cyan', 'orange', 'grey']

        for idx, var in enumerate(['rho0', 'T0', 'v02', 'v03', 'B02', 'B03']):
            equil_values = self.data.equil_data[var]
            if (equil_values == 0).all():
                continue
            self.eq_ax.plot(self.data.grid_gauss, equil_values, 'o', markersize=2, color=colors[idx], label=var)
            self.eq_ax.plot(self.data.grid_gauss, equil_values, color=colors[idx], alpha=0.3)
        self.eq_ax.set_xlabel('Gaussian grid')
        self.eq_ax.set_title('Equilibrium configuration')
        self.eq_ax.set_xlim([self.data.grid_gauss[0], self.data.grid_gauss[-1]])
        self.eq_ax.legend(loc='best')

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
                legend_item, = self.ax.plot([lb], 0, marker='p', markersize=8, color=facecolor, linestyle='none',
                                            alpha=alpha_point, label=legend_lbl)
            legend_items.append(legend_item)

        def draw_quasimodes():
            rho1 = self.data.equil_data['rho0'][0]
            rho2 = self.data.equil_data['rho0'][-1]
            a = 0.5 * self.data.params['r0']  # distance from x0 to l (linear region inbetween)
            k = np.sqrt(self.data.params['k2']**2 + self.data.params['k3']**2)
            wR_pos = np.sqrt((rho1 * wA_pos[0]**2 + rho2 * wA_pos[-1]**2) / (rho1 + rho2))
            wR_neg = np.sqrt((rho1 * wA_neg[0]**2 + rho2 * wA_neg[-1]**2) / (rho1 + rho2))
            # minus sign here, since imaginary part of the quasimode goes like exp(-wI*t)
            wI_pos = -(1.0 / 8) * np.pi * k * a * (wA_pos[-1]**2 - wA_pos[0]**2) / wR_pos
            wI_neg = -(1.0 / 8) * np.pi * k * a * (wA_neg[-1]**2 - wA_neg[0]**2) / wR_neg
            quasimodes = np.array([wR_pos + wI_pos * 1j, -wR_neg + wI_neg * 1j])
            legend_item, = self.ax.plot(np.real(quasimodes), np.imag(quasimodes), 'rs', fillstyle='none',
                                        markersize=10, linestyle='none', label='Predicted quasimodes')
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
        # if we're plotting quasimodes, annotate those as well
        if self.data.current_eq == 'resonant_absorption':
            draw_quasimodes()
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
            l_item.set_picker(10)
            region.set_visible(False)   # make continuum regions invisible by default
            regions[l_item] = region
        self.fig.canvas.mpl_connect('pick_event', on_legend_pick)
        self._plot_continua(wS_pos, wS_neg, wA_pos, wA_neg, wTH)

    def _plot_continua(self, wS_pos, wS_neg, wA_pos, wA_neg, wTH):
        self.co_ax[0].plot(self.data.grid_gauss, wS_pos, color='red', linestyle='solid',
                           label=r'$\Omega_S^-$ slow continuum')
        self.co_ax[0].plot(self.data.grid_gauss, wS_neg, color='red', linestyle='dashed',
                           label=r'$\Omega_S^+$ slow continuum')
        self.co_ax[1].plot(self.data.grid_gauss, wA_pos, color='blue', linestyle='solid',
                           label=r'$\Omega_A^-$ Alfvén continuum')
        self.co_ax[1].plot(self.data.grid_gauss, wA_neg, color='blue', linestyle='dashed',
                           label=r'$\Omega_A^+$ Alfvén continuum')
        self.co_ax[2].plot(self.data.grid_gauss, wTH, color='green', linestyle='solid',
                           label=r'Im$(\Omega_{th})$ thermal continuum')
        self.co_ax[2].set_xlabel('Gaussian grid')
        for ax_panel in self.co_ax:
            ax_panel.legend(loc='best')
            ax_panel.set_xlim([self.data.grid_gauss[0], self.data.grid_gauss[-1]])
        self.co_fig.subplots_adjust(hspace=0.05)
