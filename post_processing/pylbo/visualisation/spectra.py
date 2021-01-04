from pylbo.visualisation.figure_manager import FigureWindow


class SingleSpectrumPlot(FigureWindow):
    def __init__(self, dataset, figsize, custom_figure, **kwargs):
        super().__init__(
            figure_type="single-spectrum",
            figsize=figsize,
            custom_figure=custom_figure
        )
