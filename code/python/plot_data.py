import numpy as np
import os, sys
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import read_data

def plot_matrix(fig, ax, ax_idx, matrix, title=None):
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

            if abs(m - 0.0) < 1.0e-12:
                continue
            else:
                elements.append(m)
                row_idxs.append(i+1)    # for visual and fortran indexing
                col_idxs.append(j+1)

    im = ax[ax_idx].scatter(col_idxs, row_idxs, s=16, c=elements, cmap='jet',
                            norm=mpl.colors.LogNorm())

    ax[ax_idx].invert_yaxis()

    if title is not None:
        ax[ax_idx].set_title(title)
    ax[ax_idx].set_aspect("equal")

    divider = make_axes_locatable(ax[ax_idx])
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(im, cax=cax)


def plot_spectrum(fig, ax, omegas):
    ax.plot(np.real(omegas), np.imag(omegas), '.b', alpha=0.8)
    ax.set_xlabel(r"Re($\omega$)")
    ax.set_ylabel(r"Im($\omega$)")
    ax.set_title("Spectrum")


if __name__ == '__main__':
    config_dict = read_data.read_config_file("output/config.txt")
    gridpts = config_dict["Gridpoints"]
    matrix_gridpts = config_dict["Matrix gridpoints"]
    eigenf_gridpts = config_dict["Eigenfunction gridpoints"]

    filename_grid = "output/grid.dat"
    filename_matA = "output/matrix_A.dat"
    filename_matB = "output/matrix_B.dat"
    filename_spec = "output/eigenvalues.dat"

    eigenf_rho    = "output/eigenfunctions/rho_eigenfunction.dat"

    grid = read_data.read_stream_data(filename_grid, content_type=np.float64,
                    rows=1, cols=eigenf_gridpts)
    omegas   = read_data.read_stream_data(filename_spec, content_type=np.complex,
                                          rows=1, cols=matrix_gridpts)[0]

    matrix_B = read_data.read_stream_data(filename_matB, content_type=np.float64,
                    rows=matrix_gridpts, cols=matrix_gridpts)
    matrix_A = read_data.read_stream_data(filename_matA, content_type=np.complex,
                    rows=matrix_gridpts, cols=matrix_gridpts)

    fig, ax = plt.subplots(1, figsize=(12, 8))
    plot_spectrum(fig, ax, omegas)
    fig.tight_layout()
    fig, ax = plt.subplots(1, 2, figsize=(12, 8))
    plot_matrix(fig, ax, 1, matrix_B, title="Matrix B")
    plot_matrix(fig, ax, 0, matrix_A, title="Matrix A")
    fig.tight_layout()

    plt.show()
