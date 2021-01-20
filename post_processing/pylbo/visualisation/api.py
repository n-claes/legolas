from pylbo.data_containers import LegolasDataSet, LegolasDataSeries
from pylbo.visualisation.spectra import SingleSpectrumPlot
from pylbo.visualisation.spectra import MultiSpectrumPlot
from pylbo.utilities.logger import pylboLogger

forbidden_args = ["linestyle", "linewidth", "lw"]


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
        If the argument `dataset` is not an instance of
        :class:~pylbo.data_containers.LegolasDataSet.

    Returns
    -------
    p : `~pylbo.visualisation.figure_manager.spectra.SingleSpectrumPlot` instance
        The spectrum instance which can be used further to add continua,
        eigenfunctions, etc.
    """
    if not isinstance(dataset, LegolasDataSet):
        raise TypeError("plot_spectrum needs a single dataset, not a series.")
    for arg in forbidden_args:
        if kwargs.pop(arg, None) is not None:
            pylboLogger.warning(f"plot_spectrum does not accept the '{arg}' argument.")
    p = SingleSpectrumPlot(dataset, figsize, custom_figure, **kwargs)
    return p


def plot_spectrum_multi(
    dataseries,
    xdata,
    use_squared_omega=False,
    use_real_parts=True,
    figsize=None,
    custom_figure=None,
    **kwargs,
):
    """
    Plots the spectrum of a dataset series.

    Parameters
    ----------
    dataseries : `~pylbo.data_containers.LegolasDataSeries` instance
        The dataseries that should be used.
    xdata : str, list, np.ndarray
        Data to use for the horizontal axis. This can either be a key from the
        parameters dictionary, or a list/numpy array containing actual data.
    use_squared_omega : bool
        If `True`, this will square the eigenvalues when they are plotted on the
        vertical axis. If `False` (default), either the real or imaginary part of the
        eigenvalues will be plotted depending on the value of `use_real_parts`.
    use_real_parts : bool
        If `True` (default), this will plot the real part of the eigenvalues on the
        vertical axis. If `False` the imaginary part will be used.
    figsize : tuple
        Optional figure size like the usual matplotlib (x, x) size.
    custom_figure : tuple
        Optional, in the form (fig, ax). If supplied no new figure will be created
        but this one will be used instead. `fig` refers to the matplotlib figure and
        `ax` to a (single) axes instance, meaning that you can pass a subplot as well.

    Raises
    ------
    TypeError
        If the argument `dataseries` is not an instance of
        :class:`pylbo.data_containers.LegolasDataSeries`.

    Returns
    -------
    p : `~pylbo.visualisation.figure_manager.spectra.MultiSpectrumPlot` instance
        The spectrum instance which can be used further to add continua,
        eigenfunctions, etc.
    """
    if isinstance(dataseries, list):
        for ds in dataseries:
            if not isinstance(ds, LegolasDataSet):
                raise TypeError("invalid dataset passed to plot_spectrum_multi.")
    elif not isinstance(dataseries, LegolasDataSeries):
        raise TypeError("plot_spectrum_multi needs a dataseries, not a single dataset.")
    for arg in forbidden_args:
        if kwargs.pop(arg, None) is not None:
            pylboLogger.warning(
                f"plot_spectrum_multi does not accept the '{arg}' argument."
            )
    p = MultiSpectrumPlot(
        dataseries,
        xdata,
        use_squared_omega,
        use_real_parts,
        figsize,
        custom_figure,
        **kwargs,
    )
    return p
