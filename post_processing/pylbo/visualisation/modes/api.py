from typing import Union

import numpy as np
from pylbo.data_containers import LegolasDataSet, ensure_dataset
from pylbo.utilities.toolbox import transform_to_numpy
from pylbo.visualisation.modes.mode_data import ModeVisualisationData
from pylbo.visualisation.modes.spatial_2d import SpatialCartesianPlot2D
from pylbo.visualisation.modes.temporal_1d import TemporalEvolutionPlot1D


def plot_1d_temporal_evolution(
    ds: LegolasDataSet,
    omega: complex,
    ef_name: str,
    u2: float,
    u3: float,
    time: Union[list, np.ndarray],
    figsize: tuple[int, int] = None,
    add_background: bool = False,
    use_real_part: bool = True,
    complex_factor: complex = None,
    **kwargs,
) -> TemporalEvolutionPlot1D:
    """
    Plot the temporal evolution of a 1D eigenmode solution, i.e.
    .. math::
       f(u_1, u_2, u_3, t) = f_1(u_1) \\exp\\left(ik_2u_2 + ik_3u_3 - i\\omega t\\right)
    for a particular set of coordinates :math:`(u_2, u_3)` over a time interval.

    Parameters
    ----------
    ds : LegolasDataSet
        The data set containing the eigenfunction.
    omega : complex
        The (approximate) eigenvalue of the mode to visualise.
    ef_name : str
        The name of the eigenfunction to visualise.
    u2 : float
        The y or :math:`\\theta` coordinate of the eigenmode solution.
    u3 : float
        The z coordinate of the eigenmode solution.
    time : list or np.ndarray
        The time interval to visualise.
    figsize : tuple[int, int]
        The size of the figure.
    add_background : bool
        Whether to add the equilibrium background to the plot.
    use_real_part : bool
        Whether to use the real part of the eigenmode solution.
    complex_factor : complex
        A complex factor to multiply the eigenmode solution with.
    kwargs : dict
        Additional keyword arguments to pass to :meth:`~matplotlib.pyplot.imshow`.

    Returns
    -------
    TemporalEvolutionPlot1D
        The plot object.
    """
    ensure_dataset(ds)
    data = ModeVisualisationData(
        ds, omega, ef_name, use_real_part, complex_factor, add_background
    )
    time = transform_to_numpy(time)
    p = TemporalEvolutionPlot1D(data, u2, u3, time, figsize, **kwargs)
    p.draw()
    return p


def plot_2d_slice(
    ds: LegolasDataSet,
    omega: complex,
    ef_name: str,
    u2: Union[float, np.ndarray],
    u3: Union[float, np.ndarray],
    time: float,
    slicing_axis: str = "z",
    figsize: tuple[int, int] = None,
    add_background: bool = False,
    use_real_part: bool = True,
    complex_factor: complex = None,
    **kwargs,
) -> SpatialCartesianPlot2D:
    """
    Plot a 2D spatial view of the eigenmode solution, i.e.
    .. math::
       f(u_1, u_2, u_3, t) = f_1(u_1) \\exp\\left(ik_2u_2 + ik_3u_3 - i\\omega t\\right)
    for a particular set of coordinates. If `slicing_axis = 'z'` then a 2D view is
    created for a given slicing point along the 'z'-axis (hence a 'xy'-view), for
    `slicing_axis = 'y'` a 'xz'-view will be created. The given spatial coordinates
    `u2` and `u3` must be consistent with the slicing axis.

    Parameters
    ----------
    ds : LegolasDataSet
        The data set containing the eigenfunction.
    omega : complex
        The (approximate) eigenvalue of the mode to visualise.
    ef_name : str
        The name of the eigenfunction to visualise.
    u2 : float
        The y or :math:`\\theta` coordinate of the eigenmode solution.
    u3 : float
        The z coordinate of the eigenmode solution.
    time : float
        The time at which to visualise the eigenmode solution.
    slicing_axis : str
        The axis to slice the 2D view along, either 'z' or 'y' with default 'z'.
    figsize : tuple[int, int]
        The size of the figure.
    add_background : bool
        Whether to add the equilibrium background to the plot.
    use_real_part : bool
        Whether to use the real part of the eigenmode solution.
    complex_factor : complex
        A complex factor to multiply the eigenmode solution with.
    kwargs : dict
        Additional keyword arguments to pass to :meth:`~matplotlib.pyplot.imshow`.

    Returns
    -------
    SpatialCartesianPlot2D
        The plot object.
    """
    ensure_dataset(ds)
    data = ModeVisualisationData(
        ds, omega, ef_name, use_real_part, complex_factor, add_background
    )
    if ds.geometry == "Cartesian":
        p = SpatialCartesianPlot2D(data, u2, u3, time, slicing_axis, figsize, **kwargs)
    p.draw()
    return p
