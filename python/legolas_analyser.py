from tkinter import filedialog
from argparse import ArgumentParser
from utilities.read_data import read_config_files
from visualisations.single_spectrum import SingleSpectrum
from visualisations.multirun_spectrum import MultiSpectrum

import tkinter as tk
try:
    import matplotlib
    matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as plt
except ImportError:
    import matplotlib
    import matplotlib.pyplot as plt
import sys


def _main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--config_file', dest='config_file')
    args = parser.parse_args()

    # Check if configuration file is supplied (eg. directly from within Legolas)
    # If so, use it, and plot a single spectrum
    config_file = args.config_file

    selected_files = []
    if config_file is None:
        root = tk.Tk()
        root.withdraw()
        selected_files = list(filedialog.askopenfilenames(parent=root, title='Select .nml file(s)'))
        root.destroy()

        for file in selected_files:
            if '.nml' not in file:
                print('File with a different extension encountered! Select .nml file(s)')
                sys.exit(1)
    else:
        selected_files.append(config_file)

    namelists = read_config_files(selected_files)

    if len(namelists) == 1:
        ss = SingleSpectrum(namelists)
        ss.plot()
    else:
        ms = MultiSpectrum(namelists)
        ms.plot()


if __name__ == '__main__':
    _main()

    plt.show()
