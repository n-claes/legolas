from pylbo.visualisation.spectra import SingleSpectrumPlot


def plot_spectrum(dataset, **kwargs):
    p = SingleSpectrumPlot(dataset, **kwargs)
    return p
