import matplotlib.pyplot as plt


class FigureContainer(dict):
    def __init__(self):
        super().__init__()
        self.stack_is_enabled = True

    def add(self, figure):
        if figure.figure_id in self.figure_id_list:
            raise ValueError(
                f"id = '{figure.figure_id}' already in existing list: "
                f"{self.figure_id_list}"
            )
        self.update({figure.figure_id: figure})

    def pop(self, figure_id):
        self._validate_figure_id(figure_id)
        plt.close(figure_id)
        super().pop(figure_id)

    def _validate_figure_id(self, figure_id):
        if figure_id not in self.figure_id_list:
            raise ValueError(
                f"id='{figure_id}' not in existing list {self.figure_id_list}"
            )

    def disable_stack(self):
        """
        This disables the figurestack by closing the figures: this prevents them
        from showing when calling plt.show(). Later on we reconstruct the
        figuremanagers.
        """
        if self.stack_is_enabled:
            return
        [figure.disable() for figure in self.figure_list]
        self.stack_is_enabled = False

    def enable_stack(self):
        """
        This enables the figurestack so figures are again accessible through
        the plt.show() methods. When we lock the stack those figures are closed,
        which essentially destroys the figuremanager but keeps the figure references
        alive. Here we simply reconstruct the figuremanager for every disabled figure
        in the stack, allowing for multiple calls to show(), save(), etc.
        """
        if not self.stack_is_enabled:
            return
        [figure.enable() for figure in self.figure_list]
        self.stack_is_enabled = True

    @property
    def number_of_figures(self):
        return len(self)

    @property
    def figure_id_list(self):
        return list(self.keys())

    @property
    def figure_list(self):
        return list(self.values())

    @property
    def is_empty(self):
        return len(self) == 0


class FigureWindow:
    figure_stack = FigureContainer()

    def __init__(self, figure_type, figsize=None, custom_figure=None):
        if custom_figure is not None:
            plt.close(custom_figure[0])
            self.fig, self.ax = custom_figure
            self.figsize = tuple(self.fig.get_size_inches())
            if self.fig.get_label() == "":
                self.figure_id = self._generate_figure_id(figure_type)
            else:
                self.figure_id = self._generate_figure_id(self.fig.get_label())
        else:
            self.figsize = figsize
            if self.figsize is None:
                self.figsize = (12, 8)
            self.figure_id = self._generate_figure_id(figure_type)
            self.fig = plt.figure(self.figure_id, figsize=self.figsize)
            self.ax = self.fig.add_subplot(111)
        self.enabled = False
        self.__class__.figure_stack.add(self)

    @classmethod
    def _generate_figure_id(cls, figure_type):
        # count occurences of this type of id in the list
        occurences = sum(
            figure_type in fig_id for fig_id in cls.figure_stack.figure_id_list
        )
        suffix = 1 + occurences
        return f"{figure_type}-{suffix}"

    def disable(self):
        if not self.enabled:
            return
        self.enabled = False
        plt.close(self.figure_id)

    def enable(self):
        if self.enabled:
            return
        self.enabled = False
        manager = plt.figure(self.figure_id, self.figsize).canvas.manager
        manager.canvas.figure = self.fig
        self.fig.set_canvas(manager.canvas)

    def set_x_scaling(self):
        pass

    def set_y_scaling(self):
        pass

    def add_continua(self):
        pass

    def add_eigenfunctions(self):
        pass

    def show(self):
        """
        Shows the current figure.
        """
        self.__class__.figure_stack.disable_stack()
        self.enable()
        plt.show()
        self.__class__.figure_stack.enable_stack()

    @classmethod
    def showall(cls, reconstruct=False):
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
        if cls.figure_stack.is_empty:
            raise ValueError("No active figures present, stack is empty.")
        if not cls.figure_stack.stack_is_enabled:
            raise ValueError("Figure stack is disabled.")
        plt.show()
        if reconstruct:
            cls.figure_stack.enable_stack()
        else:
            cls.figure_stack.clear()

    def save(self, filename):
        if not self.enabled:
            raise ValueError("This figure is disabled.")
        self.__class__.figure_stack.disable_stack()
        self.enable()
        self.fig.save(filename)
        self.__class__.figure_stack.enable_stack()

    def clear(self):
        self.ax.clear()
