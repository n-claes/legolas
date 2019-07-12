import numpy as np
import matplotlib.pyplot as plt
import os, sys

def plot_spectrum(fig, ax, marker, opacity, file):
    w_real = []
    w_imag = []
    for line in file:
        # split line at delimiter
        line = line.split(",")
        # strip line
        line = [l.strip() for l in line]

        w_real.append(float(line[0]))
        w_imag.append(float(line[1]))

    ax.plot(w_real, w_imag, marker, alpha=opacity)


if __name__ == '__main__':
    filename_code = "output/tests/test_homogeneous_code.txt"
    filename_test = "output/tests/test_homogeneous_test.txt"
    file_code = open(filename_code, 'r')
    file_test = open(filename_test, 'r')

    print(">> PYTHON: Test: comparing spectra...")


    fig, ax = plt.subplots(1, figsize=(16,8))

    plot_spectrum(fig, ax, '.b', 0.5, file_code)
    plot_spectrum(fig, ax, '.r', 0.3, file_test)
    ax.set_xlabel(r"Re($\omega$)")
    ax.set_ylabel(r"Im($\omega$)")

    save = "output/tests/figures/spectrum_TEST.png"
    plt.show()
    #fig.savefig(save)
    #print(">>         figure saved to {}".format(save))
