from typing import Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pylbo.utilities.logger import pylboLogger
from pylbo.visualisation.figure_window import FigureWindow
from pylbo.visualisation.modes.mode_data import ModeVisualisationData
from pylbo.visualisation.utils import add_axis_label


class ModeFigure(FigureWindow):
    """
    Main class to hold the figure, axes and colorbar for eigenmode visualisations.

    Parameters
    ----------
    figsize : tuple[int, int]
        The size of the figure.
    data : ModeVisualisationData
        The data used for eigenmode visualisations.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
        The figure.
    axes : dict[str, matplotlib.axes.Axes]
        The axes.
    cbar : matplotlib.colorbar.Colorbar
        The colorbar.
    cbar_ax : matplotlib.axes.Axes
        The axes for the colorbar.
    data : ModeVisualisationData
        Data object containing all data associated with the selected eigenmode.
    u1_data : np.ndarray
        The data for the :math:`u_1` coordinate.
    u2_data : Union[float, np.ndarray]
        The data for the :math:`u_2` coordinate.
    u3_data : Union[float, np.ndarray]
        The data for the :math:`u_3` coordinate.
    ef_data : Union[complex, np.ndarray]
        The data for the eigenfunction.
    t_data : Union[float, np.ndarray]
        The data for the time.
    omega_txt: matplotlib.text.Text
        The text for the :math:`\\omega` label.
    k2k3_txt: matplotlib.text.Text
        The text for the :math:`k_2-k_3` label.
    u2u3_txt: matplotlib.text.Text
        The text for the :math:`u_2-u_3` label.
    t_txt: matplotlib.text.Text
        The text for the time label.
    """

    def __init__(self, figsize: tuple[int, int], data: ModeVisualisationData) -> None:
        self.cbar = None
        self._cbar_hspace = 0.01

        if figsize is None:
            figsize = (14, 8)
        fig, axes = self._create_figure_layout(figsize)
        super().__init__(fig)
        self.axes = axes
        self.cbar_ax = self._create_cbar_axes()

        # Main data object
        self.data = data
        # textbox objects
        [setattr(self, f"{val}_txt", None) for val in ("omega", "k2k3", "u2u3", "t")]
        [self._ensure_attr_set(attr) for attr in ("_u1", "_u2", "_u3", "_time")]

        self.set_plot_arrays()
        [
            self._ensure_attr_set(f"{attr}_data")
            for attr in ("u1", "u2", "u3", "time", "ef")
        ]
        self._solutions = self.calculate_mode_solution(
            ef=self.ef_data, u2=self.u2_data, u3=self.u3_data, t=self.time_data
        )
        pylboLogger.info(f"eigenmode solution shape {self._solutions.shape}")

    def _ensure_attr_set(self, attr: str) -> None:
        if getattr(self, attr, None) is None:
            raise ValueError(f"attribute '{attr}' not set for {type(self)}")

    def _check_if_number(self, val: float, attr_name: str) -> float:
        if not isinstance(val, (int, float)):
            raise ValueError(f"expected a number for {attr_name} but got {type(val)}")
        return val

    def _check_if_array(self, array: np.ndarray, attr_name: str) -> np.ndarray:
        if not isinstance(array, np.ndarray):
            raise ValueError(
                f"expected a Numpy array for {attr_name} but got {type(array)}"
            )
        return array

    def set_plot_arrays(self) -> None:
        """
        Sets the arrays used for plotting. This should implement setting of
        :attr:`u1_data`, :attr:`u2_data`, :attr:`u3_data`, :attr:`t_data` and
        :attr:`ef_data`.
        """
        raise NotImplementedError()

    def calculate_mode_solution(
        self,
        ef: np.ndarray,
        u2: Union[float, np.ndarray],
        u3: Union[float, np.ndarray],
        t: Union[float, np.ndarray],
    ) -> np.ndarray:
        """
        Calculates the mode solution.

        Parameters
        ----------
        ef : np.ndarray
            The data for the eigenfunction.
        u2 : Union[float, np.ndarray]
            The data for the :math:`u_2` coordinate.
        u3 : Union[float, np.ndarray]
            The data for the :math:`u_3` coordinate.
        t : Union[float, np.ndarray]
            The data for the time.

        Returns
        -------
        np.ndarray
            The mode solution.
        """
        return self.data.get_mode_solution(ef=ef, u2=u2, u3=u3, t=t)

    @property
    def ax(self) -> Axes:
        """
        Returns
        -------
        matplotlib.axes.Axes
            Alias for the axes containing the eigenmode solution view.
        """
        return self.axes["view"]

    @property
    def solutions(self) -> np.ndarray:
        """
        Returns
        -------
        np.ndarray
            The solutions for the eigenmode
        """
        return self._solutions

    def draw(self) -> None:
        self.draw_eigenfunction()
        self.draw_solution()
        self.draw_textboxes()
        self.add_axes_labels()

    def draw_solution(self) -> None:
        raise NotImplementedError()

    def draw_textboxes(self) -> None:
        self.add_omega_txt(self.axes["eigfunc"], loc="top left", outside=True)
        self.add_u2u3_txt(self.axes["eigfunc"], loc="top right", outside=True)
        self.add_k2k3_txt(self.ax, loc="bottom left", color="white", alpha=0.5)

    def draw_eigenfunction(self) -> None:
        """Draws the eigenfunction to the figure."""
        ax = self.axes["eigfunc"]
        ef = getattr(self.data.eigenfunction, self.data.part_name)
        grid = self.data.ds.ef_grid
        ax.plot(grid, ef, lw=2)
        ax.axvline(x=0, color="grey", ls="--", lw=1)
        ax.set_xlim(np.min(grid), np.max(grid))
        ax.set_ylabel(self.data._ef_name_latex)

    def add_axes_labels(self) -> None:
        self.ax.set_xlabel(self.get_view_xlabel())
        self.ax.set_ylabel(self.get_view_ylabel())
        self.cbar.set_label(self.get_view_cbar_label())

    def _create_cbar_axes(self) -> Axes:
        """
        Returns
        -------
        matplotlib.axes.Axes
            The axes for the colorbar.
        """
        box = self.ax.get_position()
        position = (box.x0 + box.width + self._cbar_hspace, box.y0)
        dims = (0.02, box.height)
        return self.fig.add_axes([*position, *dims])

    def get_view_xlabel(self) -> str:
        return self.data.ds.u1_str

    def get_view_ylabel(self) -> str:
        return ""

    def get_view_cbar_label(self) -> str:
        return rf"{self.data._ef_name_latex}"

    def add_omega_txt(self, ax, **kwargs) -> None:
        """
        Creates a textbox on the axis with the value of the eigenfrequency.

        Parameters
        ----------
        ax : ~matplotlib.axes.Axes
            The axes to use for the textbox.
        **kwargs
            Additional keyword arguments to pass to :meth:`add_axis_label`.
        """
        self.omega_txt = add_axis_label(
            ax, rf"$\omega$ = {self.data.omega:.5f}", **kwargs
        )

    def add_k2k3_txt(self, ax, **kwargs) -> None:
        """
        Creates a textbox on the figure with the value of the k2 and k3 coordinates.

        Parameters
        ----------
        ax : ~matplotlib.axes.Axes
            The axes to use for the textbox.
        **kwargs
            Additional keyword arguments to pass to :meth:`add_axis_label`.
        """
        self.k2k3_txt = add_axis_label(
            ax,
            f"{self.data.ds.k2_str} = {self.data.k2} | "
            f"{self.data.ds.k3_str} = {self.data.k3}",
            **kwargs,
        )

    def add_u2u3_txt(self, ax, **kwargs) -> None:
        """
        Creates a textbox on the figure with the value of the :math:`u_2-u_3`
        coordinates.

        Parameters
        ----------
        ax : ~matplotlib.axes.Axes
            The axes to use for the textbox.
        **kwargs
            Additional keyword arguments to pass to :meth:`add_axis_label`.
        """
        self.u2u3_txt = add_axis_label(
            ax,
            rf"{self.data.ds.u2_str} = {self._u2} | {self.data.ds.u3_str} = {self._u3}",
            **kwargs,
        )

    def add_t_txt(self, ax, **kwargs) -> None:
        pass

    def _create_figure_layout(self, figsize: tuple[int, int]) -> tuple[Figure, dict]:
        """
        Create the figure layout for the visualisation. Two panels are created:
        the top one for the eigenfunction and the bottom one for the visualisation.

        Parameters
        ----------
        figsize : tuple[int, int]
            The size of the figure.

        Returns
        -------
        fig : ~matplotlib.figure.Figure
            The figure to use for the visualisation.
        axes : dict
            The axes to use for the visualisation.
        """
        fig = plt.figure(figsize=figsize)
        width = 0.75
        height_1 = 0.2
        height_2 = 0.5
        v_space = 0.1
        x = (1 - width) / 2
        y1 = 1 - height_1 - v_space
        y2 = v_space
        # left, bottom, width, height in Figure coordinates
        ax1 = fig.add_axes([x, y1, width, height_1])
        ax2 = fig.add_axes([x, y2, width, height_2], sharex=ax1)
        return fig, {"eigfunc": ax1, "view": ax2}
