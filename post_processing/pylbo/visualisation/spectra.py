from pylbo.visualisation.figure_manager import FigureWindow


class SingleSpectrumPlot(FigureWindow):
    def __init__(self, dataset, **kwargs):
        super().__init__(figure_type="single-spectrum")
