from data_container import LEGOLASDataContainer
from visualisations.eigenfunctions import EigenfunctionHandler
from utilities.plot_data import plot_spectrum, plot_matrix

import matplotlib.pyplot as plt
import matplotlib.patches as patches
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
        legend_handles = []
        legend_names = []

        def draw_region(array, facecolor, alpha, legend_text=None):
            lb = np.min(array)
            rb = np.max(array)

            if np.abs(lb - rb) > 1e-10:
                self.ax.axvspan(np.min(array), np.max(array), facecolor=facecolor, alpha=alpha)
            else:   # handles the case where range of frequencies is collapsed to a single point
                self.ax.plot([rb], 0, marker='X', markersize=8, color=facecolor, alpha=0.6)

            if legend_text is not None:
                rect = patches.Rectangle((0, 0), 1, 1, facecolor=facecolor, alpha=5*alpha)
                legend_handles.append(rect)
                legend_names.append(legend_text)

        k2 = self.data.params['k2']
        k3 = self.data.params['k3']
        k0 = np.sqrt(k2**2 + k3**2)

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

        # --- Cartesian case ---
        if self.data.geometry == 'Cartesian':
            # Alfven and slow genuine singularities
            wA2 = (1/rho) * (k2*B02 + k3*B03)**2
            wS2 = (gamma*p / (gamma*p + B0**2)) * wA2

            wm2 = k0**2 * (gamma*p + B0**2) / rho

        # --- cylindrical case ---
        else:
            r = self.data.grid_gauss
            F = (k2 * B02 / r) + k3 * B03
            h2 = (k2**2 / r**2) + k3**2
            # Alfven and slow genuine singularities
            wA2 = F**2 / rho
            wS2 = (gamma*p / (gamma*p + B0**2)) * F**2 / rho

            wm2 = h2 * (gamma*p + B0**2) / rho

        # slow and fast apparent singularities (Cartesian/cylindrical differ through wm2)
        ws02 = 0.5 * wm2 * (1 - np.sqrt(1 - 4 * wS2 / wm2))
        wf02 = 0.5 * wm2 * (1 + np.sqrt(1 - 4 * wS2 / wm2))

        # split plotting of continua and singularities to avoid cluttering the plot
        # plot continua on the left side
        draw_region(doppler_shift - np.sqrt(wA2), facecolor='green', alpha=0.05,
                    legend_text=r'$\Omega_A^-(x)$ Alfven continuum')
        draw_region(doppler_shift - np.sqrt(wS2), facecolor='red', alpha=0.05,
                    legend_text=r'$\Omega_S^-(x)$ slow continuum')
        # plot apparent singularities on right side
        draw_region(doppler_shift + np.sqrt(ws02), facecolor='yellow', alpha=0.05,
                    legend_text=r'$\Omega_{s0}^+(x)$ apparent slow singularity')
        draw_region(doppler_shift + np.sqrt(wf02), facecolor='blue', alpha=0.05,
                    legend_text=r'$\Omega_{f0}^+(x)$ apparent fast singularity')

        self.ax.legend(legend_handles, legend_names, loc='best')
