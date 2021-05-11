import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylbo.visualisation.figure_manager import FigureWindow


class MatrixFigure(FigureWindow):
    """Figure showing both matrices from a dataset."""

    def __init__(self, dataset, figsize, **kwargs):
        super().__init__(figure_type="matrices", figsize=figsize)
        self.dataset = dataset
        self.kwargs = kwargs

        self.ax.change_geometry(1, 2, 1)
        self.ax2 = self.fig.add_subplot(122)
        self.draw()

    def draw(self):
        """Draws the matrices."""
        super().draw()
        self._add_matrices()

    def _add_matrices(self):
        """Adds the matrices to the figure windows."""
        # matrix A
        rows, cols, vals = self.dataset.get_matrix_A()
        # take modulus of values
        vals = np.absolute(vals)
        im = self.ax.scatter(cols, rows, c=vals, s=6, cmap="plasma", norm=LogNorm())
        self.ax.set_title("Matrix A (modulus)")
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)

        # matrix B
        rows, cols, vals = self.dataset.get_matrix_B()
        im = self.ax2.scatter(cols, rows, c=vals, s=6, cmap="plasma", norm=LogNorm())
        self.ax2.set_title("Matrix B")
        divider = make_axes_locatable(self.ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        self.fig.canvas.draw()

        for ax in (self.ax, self.ax2):
            maingrid_ticks = np.arange(0.5, self.dataset.matrix_gridpoints + 16.5, 16)
            for i in maingrid_ticks:
                ax.vlines(
                    x=i,
                    color="grey",
                    alpha=0.6,
                    ymin=0.5,
                    ymax=self.dataset.matrix_gridpoints + 0.5,
                )
                ax.hlines(
                    y=i,
                    color="grey",
                    alpha=0.6,
                    xmin=0.5,
                    xmax=self.dataset.matrix_gridpoints + 0.5,
                )
            minorgrid_ticks = np.arange(0.5, self.dataset.matrix_gridpoints + 0.5, 2)
            for i in minorgrid_ticks:
                ax.vlines(
                    x=i,
                    color="grey",
                    alpha=0.1,
                    ymin=0.5,
                    ymax=self.dataset.matrix_gridpoints + 0.5,
                )
                ax.hlines(
                    y=i,
                    color="grey",
                    alpha=0.1,
                    xmin=0.5,
                    xmax=self.dataset.matrix_gridpoints + 0.5,
                )
            visualticks = np.arange(0, self.dataset.matrix_gridpoints + 0.1, 32)
            ax.set_xticks(visualticks)
            ax.set_yticks(visualticks)
            ax.set_xlim(0, self.dataset.matrix_gridpoints + 1)
            ax.set_ylim(0, self.dataset.matrix_gridpoints + 1)
            ax.tick_params(which="both", labelsize=13)
            ax.set_aspect("equal")
            ax.invert_yaxis()
