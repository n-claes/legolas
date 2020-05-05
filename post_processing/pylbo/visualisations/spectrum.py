import matplotlib.pyplot as plt
import matplotlib.lines as mpl_lines
import matplotlib.patches as mpl_patches
import numpy as np
from .eigenfunctions import EigenfunctionHandler

class PlotSpectrum:
    def __init__(self, ds, annotate_continua=True, plot_continua=True):
        self.ds = ds
        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))
        self.annotate_continua = annotate_continua

        self._plot(self.fig, self.ax)
        self.fig.canvas.draw()
        self.fig.tight_layout()
        if plot_continua:
            self._plot_continua()

    def plot_eigenfunctions(self, merge_figs=True):
        efh = EigenfunctionHandler(self, merge_figs)
        efh.connect_figure_events()

    @staticmethod
    def show():
        plt.show()

    def _plot(self, fig, ax):
        ax.plot(self.ds.eigenvals.real, self.ds.eigenvals.imag, '.b', alpha=0.8)
        ax.axhline(y=0, linestyle='dotted', color='grey', alpha=0.3)
        ax.axvline(x=0, linestyle='dotted', color='grey', alpha=0.3)
        ax.set_xlabel(r"Re($\omega$)")
        ax.set_ylabel(r"Im($\omega$)")
        if self.annotate_continua:
            self._annotate_continua(fig, ax)

    def _annotate_continua(self, fig, ax):
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
            fig.canvas.draw()

        def draw_region(array, facecolor, legend_lbl=None, draw_horizontal=False):
            lb = np.min(array)
            rb = np.max(array)
            if np.abs(lb - rb) > 1e-10:
                if draw_horizontal:
                    legend_item = ax.axhspan(np.min(array), np.max(array), facecolor=facecolor,
                                             alpha=alpha_box, label=legend_lbl)
                else:
                    legend_item = ax.axvspan(np.min(array), np.max(array), facecolor=facecolor,
                                             alpha=alpha_box, label=legend_lbl)
            else:  # handles the case where continua are collapsed to a single point
                legend_item, = ax.plot([lb], 0, marker='p', markersize=8, color=facecolor, linestyle='none',
                                       alpha=alpha_point, label=legend_lbl)
            legend_items.append(legend_item)

        # draw regions on plot
        draw_region(self.ds.continua['wS+'], facecolor='red',
                    legend_lbl=r'$\Omega_S^+$ slow continuum')
        draw_region(self.ds.continua['wS-'], facecolor='red',
                    legend_lbl=r'$\Omega_S^-$ slow continuum')
        draw_region(self.ds.continua['wA+'], facecolor='cyan',
                    legend_lbl=r'$\Omega_A^+$ Alfvén continuum')
        draw_region(self.ds.continua['wA-'], facecolor='cyan',
                    legend_lbl=r'$\Omega_A^-$ Alfvén continuum')
        # thermal continuum is imaginary, so draw horizontally
        draw_region(self.ds.continua['wth'], facecolor='green',
                    legend_lbl=r'$\Omega_{th}$ thermal continuum', draw_horizontal=True)
        legend = ax.legend(loc='best')

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
            l_item.set_alpha(0.1)
            region.set_visible(False)  # make continuum regions invisible by default
            regions[l_item] = region
        fig.canvas.mpl_connect('pick_event', on_legend_pick)

    def _plot_continua(self):
        fig, ax = plt.subplots(3, 1, figsize=(12, 8))
        ax[0].plot(self.ds.grid_gauss, self.ds.continua['wS-'], color='red',
                   linestyle='solid', label=r'$\Omega_S^-$ slow continuum')
        ax[0].plot(self.ds.grid_gauss, self.ds.continua['wS+'], color='red',
                   linestyle='dashed', label=r'$\Omega_S^+$ slow continuum')
        ax[1].plot(self.ds.grid_gauss, self.ds.continua['wA-'], color='blue',
                   linestyle='solid', label=r'$\Omega_A^-$ Alfvén continuum')
        ax[1].plot(self.ds.grid_gauss, self.ds.continua['wA+'], color='blue',
                   linestyle='dashed', label=r'$\Omega_A^+$ Alfvén continuum')
        ax[2].plot(self.ds.grid_gauss, self.ds.continua['wth'], color='green',
                   linestyle='solid', label=r'Im$(\Omega_{th})$ thermal continuum')
        ax[2].set_xlabel('Grid')
        for ax_panel in ax:
            ax_panel.legend(loc='best')
            ax_panel.set_xlim([self.ds.grid_gauss[0], self.ds.grid_gauss[-1]])
        fig.subplots_adjust(hspace=0.05)
        fig.tight_layout()
