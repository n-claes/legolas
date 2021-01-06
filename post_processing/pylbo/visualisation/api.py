from pylbo.data_containers import LegolasDataSet
from pylbo.visualisation.spectra import SingleSpectrumPlot
from pylbo.utilities.logger import pylboLogger


def plot_spectrum(dataset, figsize=None, custom_figure=None, **kwargs):
    if not isinstance(dataset, LegolasDataSet):
        raise TypeError("plot_spectrum needs a single dataset, not a series.")
    forbidden_args = ["linestyle", "linewidth", "lw"]
    for arg in forbidden_args:
        if kwargs.pop(arg, None) is not None:
            pylboLogger.warning(f"plot_spectrum does not accept the '{arg}' argument.")
    p = SingleSpectrumPlot(dataset, figsize, custom_figure, **kwargs)
    return p
