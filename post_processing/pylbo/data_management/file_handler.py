import tkinter as tk
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

def select_files():
    root = tk.Tk()
    root.withdraw()
    root.lift()
    root.focus_set()
    files = list(filedialog.askopenfilenames(parent=root, title='Select Legolas file(s)'))
    root.destroy()
    return files