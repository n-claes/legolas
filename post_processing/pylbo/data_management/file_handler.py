import tkinter as tk
import numpy as np
from pathlib import Path
from tkinter import filedialog
from .data_container import LegolasDataContainer

def load(datfiles):
    if isinstance(datfiles, list):
        datasets = []
        for datfile in datfiles:
            ds = LegolasDataContainer(datfile)
            datasets.append(ds)
        if len(datasets) == 1:
            datasets = datasets.pop()
        return datasets
    ds = LegolasDataContainer(datfiles)
    return ds

def read_log_file(file, sort=False):
    if not ".log" in str(file):
        raise ValueError('This is not a .log file: {}'.format(file))
    filepath = Path(file).resolve()
    if not filepath.is_file():
        raise FileNotFoundError(filepath)
    eigenvalues = []
    with open(filepath, 'r') as logfile:
        for line in logfile:
            line = line.strip().split(',')
            x, y = line
            eigenvalues.append(complex(float(x), float(y)))
    eigenvalues = np.asarray(eigenvalues)
    if sort:
        eigenvalues = np.sort(eigenvalues)
    return eigenvalues

def select_files():
    root = tk.Tk()
    root.withdraw()
    root.lift()
    root.focus_set()
    files = list(filedialog.askopenfilenames(parent=root, title='Select Legolas file(s)'))
    root.destroy()
    if not files:
        exit()
    return files
