import matplotlib.pyplot as plt


class FigureContainer(dict):
    def add(self, figure):
        if figure.figure_id in self.figure_id_list:
            raise ValueError(
                f"id = '{figure.figure_id}' already in existing list: "
                f"{self.figure_id_list}"
            )
        self.update({figure.figure_id: figure})

    def pop(self, figure):
        self._validate_figure_id(figure.figure_id)
        return self.pop(figure.figure_id)

    def get(self, figure):
        self._validate_figure_id(figure.figure_id)
        return self.get(figure.figure_id)

    def remove(self, figure):
        figure = self.pop(figure.figure_id)
        plt.delaxes(figure.ax)
        plt.close(figure.fig)

    def _validate_figure_id(self, figure_id):
        if figure_id not in self.figure_id_list:
            raise ValueError(
                f"id='{figure_id}' not in existing list {self.figure_id_list}"
            )

    @property
    def number_of_figures(self):
        return len(self)

    @property
    def figure_id_list(self):
        return list(self.keys())

    @property
    def is_empty(self):
        return len(self) == 0


class FigureWindow:
    figure_stack = FigureContainer()

    def __init__(self, figure_type, figsize=None):
        self.figsize = figsize
        if self.figsize is None:
            self.figsize = (12, 8)
        self.figure_id = self._generate_figure_id(figure_type)
        self.fig = plt.figure(self.figure_id, figsize=self.figsize)
        self.ax = self.fig.add_subplot(111)
        self.ax.plot([0, 1, 2], [2, 3, 4])
        self.__class__.figure_stack.add(self)

    @classmethod
    def _generate_figure_id(cls, figure_type):
        # count occurences of this type of id in the list
        occurences = sum(
            figure_type in fig_id for fig_id in cls.figure_stack.figure_id_list
        )
        suffix = 1 + occurences
        return f"{figure_type}-{suffix}"

    def use_custom_figure(self, fig, ax):
        pass

    def set_x_scaling(self):
        pass

    def set_y_scaling(self):
        pass

    def add_continua(self):
        pass

    def add_eigenfunctions(self):
        pass

    def _reconstruct_figure_managers(self):
        """
        plt.close() and plt.show() destroy the figuremanager but the figure
        references are kept alive in the figure stack. Here we simply reconstruct
        the figuremanager for every figure in the stack, allowing for multiple
        calls to show(), save(), etc.
        """
        for fig_id, figure in self.__class__.figure_stack.items():
            manager = plt.figure(fig_id, self.figsize).canvas.manager
            manager.canvas.figure = figure.fig
            figure.fig.set_canvas(manager.canvas)

    def _close_other_figures(self):
        """
        Closes all figures besides the current instance.
        """
        for fig_id, figure in self.__class__.figure_stack.items():
            if fig_id == self.figure_id:
                continue
            plt.close(figure.fig)

    def show(self):
        """
        Shows (only) the current figure, uses a workaround since matplotlib
        shows all active figures at once. We want to be able to show the figure,
        close it and then modify it and show it again, or save it.
        """
        self._close_other_figures()
        plt.show()
        self._reconstruct_figure_managers()

    def showall(self, reconstruct=False):
        """
        Shows all active figures at once through a call to plt.show().
        Unless `reconstruct=True` this is final: all figures and figuremanagers
        will be destroyed and will not be reconstructed.
        Parameters
        ----------
        reconstruct : bool
            If `True`, reconstructs all figuremanagers, allowing for multiple
            calls to the show routines.
        """
        if self.__class__.figure_stack.is_empty:
            raise ValueError("No active figures present, stack is empty.")
        plt.show()
        if reconstruct:
            self._reconstruct_figure_managers()
        else:
            self.__class__.figure_stack.clear()

    def save(self, filename):
        self._close_other_figures()
        self.fig.save(filename)
        self._reconstruct_figure_managers()

    def clear(self):
        self.ax.clear()
