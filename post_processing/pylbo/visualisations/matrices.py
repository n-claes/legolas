import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..utilities.logger import pylboLogger


def plot_matrices(ds, force_draw=False):
    """
    Plots the matrices for a given dataset. Creates one single figure with the
    A-matrix on the left and the B-matrix on the right.

    Warnings
    --------
    This method can take a long time for large datasets due to the amount of points
    that are being plotted. This method is meant solely for inspections and should
    not be used for visualisation purposes. Visualising the matrices of large datasets
    will not be clear anyway due to the sheer amount of blocks present. By default,
    a warning message will be printed when the amount of gridpoints is higher than 20
    and plotting will be skipped.
    Supply the kwarg `force_draw=True` if you want to draw the figure anyway.

    Parameters
    ----------
    ds : ~pylbo.LegolasDataContainer
        The :class:`~pylbo.LegolasDataContainer` instance currently loaded.
    force_draw : bool
        Forces drawing the matrix figure, defaults to `False`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure on which the matrixes are drawn
    axes : numpy.ndarray(dtype=matplotlib.axes.Axes, ndim=1)
        The axes on which the matrices are drawn. Index 0 corresponds to the A-matrix,
        index 1 to the B-matrix.
    """
    if ds.gridpts > 20 and not force_draw:
        pylboLogger.warning(
            "plotting the matrices is disabled if gridpoints > 20, "
            "use force_draw=True to override."
        )
        return None, None
    pylboLogger.info("drawing matrix figure...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axl, axr = axes.flatten()

    # matrix A
    rows, cols, vals = ds.get_matrix_A()
    # take modulus of values
    vals = np.absolute(vals)
    im = axl.scatter(cols, rows, c=vals, s=6, cmap="plasma", norm=LogNorm())
    axl.set_title("Matrix A (modulus)")
    divider = make_axes_locatable(axl)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    # matrix B
    rows, cols, vals = ds.get_matrix_B()
    im = axr.scatter(cols, rows, c=vals, s=6, cmap="plasma", norm=LogNorm())
    axr.set_title("Matrix B")
    divider = make_axes_locatable(axr)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    fig.canvas.draw()

    for ax in (axl, axr):
        maingrid_ticks = np.arange(0.5, ds.matrix_gridpts + 16.5, 16)
        for i in maingrid_ticks:
            ax.vlines(
                x=i, color="grey", alpha=0.6, ymin=0.5, ymax=ds.matrix_gridpts + 0.5
            )
            ax.hlines(
                y=i, color="grey", alpha=0.6, xmin=0.5, xmax=ds.matrix_gridpts + 0.5
            )
        minorgrid_ticks = np.arange(0.5, ds.matrix_gridpts + 0.5, 2)
        for i in minorgrid_ticks:
            ax.vlines(
                x=i, color="grey", alpha=0.1, ymin=0.5, ymax=ds.matrix_gridpts + 0.5
            )
            ax.hlines(
                y=i, color="grey", alpha=0.1, xmin=0.5, xmax=ds.matrix_gridpts + 0.5
            )

        visualticks = np.arange(0, ds.matrix_gridpts + 0.1, 32)
        ax.set_xticks(visualticks)
        ax.set_yticks(visualticks)
        ax.set_xlim(0, ds.matrix_gridpts + 1)
        ax.set_ylim(0, ds.matrix_gridpts + 1)
        ax.tick_params(which="both", labelsize=13)
        ax.set_aspect("equal")
        ax.invert_yaxis()

    pylboLogger.info("done!")

    return fig, axes
