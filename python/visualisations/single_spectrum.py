from data_container import LEGOLASDataContainer
from visualisations.eigenfunctions import EigenfunctionHandler
from utilities.plot_data import plot_spectrum, plot_matrix

import matplotlib.pyplot as plt

class SingleSpectrum(LEGOLASDataContainer):
    def __init__(self, namelist):
        super().__init__(namelist)

        self.omegas = self.read_omegas()
        self.matA = None
        self.matB = None

        self.fig, self.ax = plt.subplots(1, figsize=(12, 8))


    def plot(self):
        plot_spectrum(self.fig, self.ax, self.omegas, title=self.current_eq)

        # if self.show_mats:
        #     self._plot_matrices()
        if self.show_eigenfuncs:
            efh = EigenfunctionHandler(self)
            efh.connect_interactive_funcs()

        plt.show()


    def _plot_matrices(self):
        self.matA = self.read_matA()
        self.matB = self.read_matB()

        fig, ax = plt.subplots(1, 2, figsize=(12, 8))
        plot_matrix(fig, ax, 0, self.matA, title='Matrix A')
        plot_matrix(fig, ax, 1, self.matB, title='Matrix B')
