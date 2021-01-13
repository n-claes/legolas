from pylbo.data_containers import LegolasDataSet
from pylbo.visualisation.spectra import SingleSpectrumPlot
from pylbo.utilities.logger import pylboLogger


def plot_spectrum(dataset, figsize=None, custom_figure=None, **kwargs):
    """
    Plots the spectrum of a single dataset.

    Parameters
    ----------
    dataset : `~pylbo.data_containers.LegolasDataSet` instance
        The dataset that should be used.
    figsize : tuple
        Optional figure size like the usual matplotlib (x, x) size.
    custom_figure : tuple
        Optional, in the form (fig, ax). If supplied no new figure will be created
        but this one will be used instead. `fig` refers to the matplotlib figure and
        `ax` to a (single) axes instance, meaning that you can pass a subplot as well.

    Raises
    ------
    TypeError
        If the argument `dataset` is of a wrong type.

    Returns
    -------
    p : `~pylbo.visualisation.figure_manager.spectra.SingleSpectrumPlot` instance
        The spectrum instance which can be used further to add continua,
        eigenfunctions, etc.
    """
    if not isinstance(dataset, LegolasDataSet):
        raise TypeError("plot_spectrum needs a single dataset, not a series.")
    forbidden_args = ["linestyle", "linewidth", "lw"]
    for arg in forbidden_args:
        if kwargs.pop(arg, None) is not None:
            pylboLogger.warning(f"plot_spectrum does not accept the '{arg}' argument.")
    p = SingleSpectrumPlot(dataset, figsize, custom_figure, **kwargs)
    return p
