import matplotlib.pyplot as plt
from .data_management.api import load, \
    select_files
from .visualisations.api import PlotSpectrum, \
    MultiSpectrum

def show():
    plt.show()
