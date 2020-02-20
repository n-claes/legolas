from data_container import LEGOLASDataContainer
from visualisations.eigenfunctions import EigenfunctionHandler
from utilities.plot_data import plot_spectrum, plot_matrix

import matplotlib.pyplot as plt

class SingleSpectrum(LEGOLASDataContainer):
    def __init__(self, namelist_array):
        super().__init__(namelist_array)

        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))
        self.data = self.datacontainer.pop()


    def plot(self):
        plot_spectrum(self.fig, self.ax, self.data.omegas, title=self.data.current_eq)

        # if data.fname_eq is not None:
        #     self._annotate_frequencies()
        #     return
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
