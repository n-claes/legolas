import numpy as np
import matplotlib.pyplot as plt
import sys, os

def read_file(filename, matrix_type="real"):
    file = open(filename, 'r')

    rows = []
    cols = []
    omegas = []
    for line in file:
        line = line.split(",")
        line = [l.strip() for l in line]

        # correct for Fortran indexing
        rows.append(int(line[0])-1)
        cols.append(int(line[1])-1)

        if matrix_type == "real":
            tmp = float(line[2])
        elif matrix_type == "complex":
            tmp = float(line[2]) + 1j*float(line[3])
        else:
            print("wrong matrix type passed")
            sys.exit()
        omegas.append(tmp)

    return rows, cols, omegas


def reduce_size(rows, cols, omegas):
    max_val = max(rows)

    moduli = [np.sqrt(w.real**2 + w.imag**2) for w in omegas]
    # Reduce array size for plotting
    moduli_new = []
    rows_new   = []
    cols_new   = []
    for i in range(len(moduli)):
        m = moduli[i]
        if (abs(m - 0.0) < 1.0e-12):
            continue
        else:
            moduli_new.append(m)
            rows_new.append(rows[i])
            cols_new.append(cols[i])

    return rows_new, cols_new, moduli_new


def plot_matrix(fig, ax, rows, cols, moduli, idx=0, title=None):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    max_idx = max(rows)
    # Major ticks every 16, minor ticks every 8
    major_ticks = np.arange(0, max_idx + 32, 32)
    minor_ticks = np.arange(0, max_idx + 16, 16)

    # Set ticks
    ax[idx].set_xticks(major_ticks)
    ax[idx].set_xticks(minor_ticks, minor=True)
    ax[idx].set_yticks(major_ticks)
    ax[idx].set_yticks(minor_ticks, minor=True)
    # Set grid
    ax[idx].grid(which='minor', linestyle="dotted", color="grey", alpha=0.3)
    ax[idx].grid(which='major', linestyle="dotted", color="grey", alpha=0.5)

    im = ax[idx].scatter(rows, cols, s=16, c=moduli, cmap='jet')
    ax[idx].invert_yaxis()

    if title is not None:
        ax[idx].set_title(title)
    ax[idx].set_aspect("equal")

    divider = make_axes_locatable(ax[idx])
    cax = divider.append_axes("right", size="5%", pad=0.05)

    plt.colorbar(im, cax=cax)


if __name__ == '__main__':
    filenameA = "output/matrix_A.txt"
    filenameB = "output/matrix_B.txt"

    print(">> PYTHON: Plotting matrices...")

    filenames = [filenameA, filenameB]
    titles = ["Matrix A", "Matrix B"]

    fig, ax = plt.subplots(1, 2, figsize=(16, 8))

    for filename in filenames:
        idx = filenames.index(filename)
        rows, cols, omegas = read_file(filename)
        rows, cols, moduli = reduce_size(rows, cols, omegas)
        plot_matrix(fig, ax, rows, cols, moduli, idx=idx, title=titles[idx])

    save = "output/figures/matrices.png"
    plt.show()
    # fig.savefig(save)
    # print(">>         figure saved to {}".format(save))
