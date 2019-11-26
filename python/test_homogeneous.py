import numpy as np
import matplotlib.pyplot as plt
import os, sys
import read_data
import plot_data




if __name__ == '__main__':
    tests = ["adiabatic", "gravity"]
    gridpts = 31
    matrix_gridpts = 16*gridpts

    for name in tests:
        filename_code   = "output/tests/homo_{}_test_CODE.dat".format(name)
        filename_theory = "output/tests/homo_{}_test_THEORY.dat".format(name)

        w_code   = read_data.read_stream_data(filename_code,
                                content_type=np.complex, rows=1,
                                cols=matrix_gridpts)
        w_theory = read_data.read_stream_data(filename_theory,
                                content_type=np.complex, rows=1,
                                cols=matrix_gridpts)


        fig, ax = plt.subplots(1, figsize=(16,8))
        title = "{} (test)".format(name)
        plot_data.plot_spectrum(fig, ax, w_code, marker='.b', alpha=0.8,
                                title=title)
        plot_data.plot_spectrum(fig, ax, w_theory, marker='rx', alpha=0.3)

    plt.show()
