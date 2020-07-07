import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_matrices(ds):
    """
    Plots the matrices for a given dataset. Creates one single figure with the
    A-matrix on the left and the B-matrix on the right.

    Warnings
    --------
    This method can take a long time for large datasets due to the amount of points
    that are being plotted. This method is meant solely for inspections and should
    not be used for visualisation purposes. Visualising the matrices of large datasets
    will not be clear anyway due to the sheer amount of blocks present.

    Parameters
    ----------
    ds : ~pylbo.LegolasDataContainer
        The :class:`~pylbo.LegolasDataContainer` instance currently loaded.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure on which the matrixes are drawn
    axes : numpy.ndarray(dtype=matplotlib.axes.Axes, ndim=1)
        The axes on which the matrices are drawn. Index 0 corresponds to the A-matrix,
        index 1 to the B-matrix.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # Major ticks every 32, minor ticks every 16
    major_ticks = np.arange(0, ds.matrix_gridpts + 32, 32)
    minor_ticks = np.arange(0, ds.matrix_gridpts + 16, 16)

    for ax in axes:
        # Set ticks
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)
        # Set grid
        ax.grid(which='minor', linestyle="dotted", color="grey", alpha=0.3)
        ax.grid(which='major', linestyle="dotted", color="black", alpha=0.5)

    # matrix A
    rows, cols, vals = ds.get_matrix_A()
    # take modulus of values
    vals = np.absolute(vals)
    im = axes[0].scatter(cols, rows, c=vals, s=20, cmap='plasma', norm=LogNorm())
    axes[0].set_aspect('equal')
    axes[0].set_title('Matrix A (modulus)')
    axes[0].invert_yaxis()
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    # matrix B
    rows, cols, vals = ds.get_matrix_B()
    im = axes[1].scatter(cols, rows, c=vals, s=20, cmap='plasma', norm=LogNorm())
    axes[1].set_aspect('equal')
    axes[1].set_title('Matrix B')
    axes[1].invert_yaxis()
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    fig.canvas.draw()

    return fig, axes
