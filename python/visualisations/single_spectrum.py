from data_container import LEGOLASDataContainer
from visualisations.eigenfunctions import EigenfunctionHandler
from utilities.plot_data import plot_spectrum, plot_matrix

import matplotlib.pyplot as plt

class SingleSpectrum(LEGOLASDataContainer):
    def __init__(self, namelist_array):
        super().__init__(namelist_array)

        self[0].fig, self[0].ax = plt.subplots(1, figsize=(12, 8))


    def plot(self):
        plot_spectrum(self[0].fig, self[0].ax, self[0].omegas, title=self[0].current_eq)

        if self[0].show_mats:
            self._plot_matrices()
        if self[0].show_equil:
            self._plot_equilibria()
        if self[0].show_eigenfuncs:
            efh = EigenfunctionHandler(self[0])
            efh.connect_interactive_funcs()


    def _plot_matrices(self):
        fig, ax = plt.subplots(1, 2, figsize=(12, 8))
        plot_matrix(fig, ax, 0, self[0].matA, title='Matrix A')
        plot_matrix(fig, ax, 1, self[0].matB, title='Matrix B')


    def _plot_equilibria(self):
        fig, ax = plt.subplots(1, figsize=(12, 8))

        for var in ['rho0', 'T0', 'v02', 'v03', 'B02', 'B03']:
            ax.plot(self[0].grid_gauss, self[0].equil_data[var], label=var)
        ax.set_xlabel('Gaussian grid')
        ax.set_title('Equilibrium configuration')
        ax.set_xlim([self[0].grid_gauss[0], self[0].grid_gauss[-1]])
        ax.legend(loc='best')
