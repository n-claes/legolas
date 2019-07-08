import numpy as np
import matplotlib.pyplot as plt
import os, sys

if __name__ == '__main__':
    filename = "eigenvalues.txt"
    file = open(filename, 'r')

    w_real = []
    w_imag = []
    for line in file:
        line = line.strip()
        line = line.split()

        w_real.append(float(line[1]))
        w_imag.append(float(line[3]))

    fig, ax = plt.subplots(1, figsize=(16,8))
    ax.plot(w_real, w_imag, '.b')
    ax.set_xlabel(r"Re($\omega$)")
    ax.set_ylabel(r"Im($\omega$)")

    plt.show()
