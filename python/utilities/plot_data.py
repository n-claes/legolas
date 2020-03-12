import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_matrix(fig, ax, ax_idx, matrix, title=None, log=True):
    matrix_gridpts = matrix.shape[0]

    # Major ticks every 16, minor ticks every 8
    major_ticks = np.arange(0, matrix_gridpts + 32, 32)
    minor_ticks = np.arange(0, matrix_gridpts + 16, 16)

    # Set ticks
    ax[ax_idx].set_xticks(major_ticks)
    ax[ax_idx].set_xticks(minor_ticks, minor=True)
    ax[ax_idx].set_yticks(major_ticks)
    ax[ax_idx].set_yticks(minor_ticks, minor=True)
    # Set grid
    ax[ax_idx].grid(which='minor', linestyle="dotted", color="grey", alpha=0.3)
    ax[ax_idx].grid(which='major', linestyle="dotted", color="grey", alpha=0.5)

    elements = []
    row_idxs = []
    col_idxs = []

    for i in range(0, matrix_gridpts):
        for j in range(0, matrix_gridpts):
            m = np.sqrt(matrix[i, j].real**2 + matrix[i, j].imag**2)

            if abs(m) >= 1.0e-12:
                elements.append(m)
                row_idxs.append(i+1)    # for visual and fortran indexing
                col_idxs.append(j+1)

    if log:
        im = ax[ax_idx].scatter(col_idxs, row_idxs, s=16, c=elements, cmap='jet', norm=mpl.colors.LogNorm())
    else:
        im = ax[ax_idx].scatter(col_idxs, row_idxs, s=16, c=elements, cmap='jet')

    ax[ax_idx].invert_yaxis()

    if title is not None:
        ax[ax_idx].set_title(title)
    ax[ax_idx].set_aspect("equal")

    divider = make_axes_locatable(ax[ax_idx])
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(im, cax=cax)
    fig.tight_layout()


def plot_spectrum(fig, ax, omegas, marker=None, alpha=0.8, title=None):
    if marker is None:
        marker = ".b"

    # Do not plot ZBIG values but keep them in original array
    omegas_plot = omegas[np.where(np.imag(omegas) < 1.0e15)]

    ax.plot(np.real(omegas_plot), np.imag(omegas_plot), marker, alpha=alpha)
    ax.axhline(y=0, linestyle='dotted', color='grey', alpha=0.3)
    ax.axvline(x=0, linestyle='dotted', color='grey', alpha=0.3)
    ax.set_xlabel(r"Re($\omega$)")
    ax.set_ylabel(r"Im($\omega$)")
    if title is not None:
        ax.set_title("Spectrum: {}".format(title))
    fig.tight_layout()


def plot_eigenfunctions(fig, ax, omegas, grid, eigenfunctions, var, w_idx, real=True):
    eigenf = eigenfunctions[var]
    if isinstance(w_idx, list):
        for w in w_idx:
            lab = r"$\omega${} = {:.8f}".format(w, omegas[w])
            if real:
                ax.plot(grid, np.real(eigenf[:, w]), label=lab)
            else:
                ax.plot(grid, np.imag(eigenf[:, w]), label=lab)
    else:
        lab = r"$\omega${} = {}".format(w_idx, omegas[w_idx])
        if real:
            ax.plot(grid, np.real(eigenf[:, w_idx]), label=lab)
        else:
            ax.plot(grid, np.imag(eigenf[:, w_idx]), label=lab)

    ax.axhline(y=0, linestyle='dotted', color='grey')
    ax.axvline(x=0, linestyle='dotted', color='grey')
    if real:
        ax.set_title(var + " (real part)")
    else:
        ax.set_title(var + " (imaginary part)")
    ax.legend(loc='best')
    fig.tight_layout()
