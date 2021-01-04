from pylbo.visualisation.spectra import SingleSpectrumPlot


def plot_spectrum(dataset, figsize=None, custom_figure=None, **kwargs):
    p = SingleSpectrumPlot(dataset, figsize, custom_figure, **kwargs)
    return p
