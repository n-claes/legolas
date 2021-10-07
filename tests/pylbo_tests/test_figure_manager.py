import pytest
import matplotlib.pyplot as plt
from pylbo.visualisation.figure_manager import FigureWindow
from pylbo.utilities.toolbox import get_axis_geometry


def test_create_window():
    p = FigureWindow("test")
    assert p.figure_id == "test-1"
    assert get_axis_geometry(p.ax) == (1, 1, 0)


def test_add_subplot_invalid():
    p = FigureWindow("test")
    with pytest.raises(ValueError):
        p._add_subplot_axes(p.ax, loc="unknown")


# ========== TESTS FOR RIGHT SUBPLOT ADDITION ==========
def test_add_subplot_singleplot_right():
    p = FigureWindow("test")
    new_ax = p._add_subplot_axes(p.ax, loc="right")
    assert get_axis_geometry(p.ax) == (1, 2, 0)
    assert get_axis_geometry(new_ax) == (1, 2, 1)
    assert len(p.fig.get_axes()) == 2


def test_add_subplot_multiple_columns_right_end():
    fig, (ax0, ax1) = plt.subplots(1, 2)
    p = FigureWindow("test", custom_figure=(fig, ax1))
    new_ax = p._add_subplot_axes(p.ax, loc="right")
    assert get_axis_geometry(ax0) == (1, 2, 0)
    assert get_axis_geometry(ax1) == (1, 2, 0)
    assert get_axis_geometry(new_ax) == (1, 2, 1)
    assert len(p.fig.get_axes()) == 3


def test_add_subplot_multiple_columns_right_middle():
    fig, (ax0, ax2) = plt.subplots(1, 2)
    p = FigureWindow("test", custom_figure=(fig, ax0))
    new_ax = p._add_subplot_axes(p.ax, loc="right")
    assert get_axis_geometry(ax0) == (1, 2, 0)
    assert get_axis_geometry(new_ax) == (1, 2, 1)
    assert get_axis_geometry(ax2) == (1, 2, 1)
    assert len(p.fig.get_axes()) == 3


def test_add_subplot_multiple_rows_right_top():
    fig, (ax0, ax1) = plt.subplots(2, 1)
    p = FigureWindow("test", custom_figure=(fig, ax0))
    new_ax = p._add_subplot_axes(p.ax, loc="right")
    assert get_axis_geometry(ax0) == (1, 2, 0)
    assert get_axis_geometry(new_ax) == (1, 2, 1)
    assert get_axis_geometry(ax1) == (2, 1, 1)
    assert len(fig.get_axes()) == 3


def test_add_subplot_multiple_rows_right_middle():
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1)
    p = FigureWindow("test", custom_figure=(fig, ax1))
    new_ax = p._add_subplot_axes(p.ax, loc="right")
    assert get_axis_geometry(ax0) == (3, 1, 0)
    assert get_axis_geometry(ax1) == (1, 2, 0)
    assert get_axis_geometry(new_ax) == (1, 2, 1)
    assert get_axis_geometry(ax2) == (3, 1, 2)
    assert len(fig.get_axes()) == 4


def test_add_subplot_matrix_right():
    fig, axes = plt.subplots(3, 3)
    p = FigureWindow("test", custom_figure=(fig, axes[1, 0]))
    new_ax = p._add_subplot_axes(p.ax, loc="right")
    assert get_axis_geometry(axes[0, 0]) == (3, 3, 0)
    assert get_axis_geometry(axes[1, 1]) == (3, 3, 4)
    assert get_axis_geometry(p.ax) == (1, 2, 0)
    assert get_axis_geometry(new_ax) == (1, 2, 1)
    assert len(fig.get_axes()) == 10


# ========== TESTS FOR LEFT SUBPLOT ADDITION ==========
def test_add_subplot_singleplot_left():
    p = FigureWindow("test")
    new_ax = p._add_subplot_axes(p.ax, loc="left")
    assert get_axis_geometry(p.ax) == (1, 2, 1)
    assert get_axis_geometry(new_ax) == (1, 2, 0)
    assert len(p.fig.get_axes()) == 2


def test_add_subplot_multiple_columns_left_end():
    fig, (ax0, ax1) = plt.subplots(1, 2)
    p = FigureWindow("test", custom_figure=(fig, ax1))
    new_ax = p._add_subplot_axes(p.ax, loc="left")
    assert get_axis_geometry(ax0) == (1, 2, 0)
    assert get_axis_geometry(ax1) == (1, 2, 1)
    assert get_axis_geometry(new_ax) == (1, 2, 0)
    assert len(p.fig.get_axes()) == 3


def test_add_subplot_multiple_columns_left_front():
    fig, (ax0, ax2) = plt.subplots(1, 2)
    p = FigureWindow("test", custom_figure=(fig, ax0))
    new_ax = p._add_subplot_axes(p.ax, loc="left")
    assert get_axis_geometry(ax0) == (1, 2, 1)
    assert get_axis_geometry(new_ax) == (1, 2, 0)
    assert get_axis_geometry(ax2) == (1, 2, 1)
    assert len(p.fig.get_axes()) == 3


def test_add_subplot_multiple_rows_left_top():
    fig, (ax0, ax1) = plt.subplots(2, 1)
    p = FigureWindow("test", custom_figure=(fig, ax0))
    new_ax = p._add_subplot_axes(p.ax, loc="left")
    assert get_axis_geometry(ax0) == (1, 2, 1)
    assert get_axis_geometry(new_ax) == (1, 2, 0)
    assert get_axis_geometry(ax1) == (2, 1, 1)
    assert len(fig.get_axes()) == 3


def test_add_subplot_multiple_rows_left_middle():
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1)
    p = FigureWindow("test", custom_figure=(fig, ax1))
    new_ax = p._add_subplot_axes(p.ax, loc="left")
    assert get_axis_geometry(ax0) == (3, 1, 0)
    assert get_axis_geometry(ax1) == (1, 2, 1)
    assert get_axis_geometry(new_ax) == (1, 2, 0)
    assert get_axis_geometry(ax2) == (3, 1, 2)
    assert len(fig.get_axes()) == 4


def test_add_subplot_matrix_left():
    fig, axes = plt.subplots(3, 3)
    p = FigureWindow("test", custom_figure=(fig, axes[0, 2]))
    new_ax = p._add_subplot_axes(p.ax, loc="left")
    assert get_axis_geometry(axes[0, 1]) == (3, 3, 1)
    assert get_axis_geometry(axes[1, 2]) == (3, 3, 5)
    assert get_axis_geometry(p.ax) == (1, 2, 1)
    assert get_axis_geometry(new_ax) == (1, 2, 0)
    assert len(fig.get_axes()) == 10


# ========== TESTS FOR TOP SUBPLOT ADDITION ==========
def test_add_subplot_singleplot_top():
    p = FigureWindow("test")
    new_ax = p._add_subplot_axes(p.ax, loc="top")
    assert get_axis_geometry(p.ax) == (2, 1, 1)
    assert get_axis_geometry(new_ax) == (2, 1, 0)
    assert len(p.fig.get_axes()) == 2


def test_add_subplot_multiple_columns_top_end():
    fig, (ax0, ax1) = plt.subplots(1, 2)
    p = FigureWindow("test", custom_figure=(fig, ax1))
    new_ax = p._add_subplot_axes(p.ax, loc="top")
    assert get_axis_geometry(ax0) == (1, 2, 0)
    assert get_axis_geometry(ax1) == (2, 1, 1)
    assert get_axis_geometry(new_ax) == (2, 1, 0)
    assert len(p.fig.get_axes()) == 3


def test_add_subplot_multiple_columns_top_front():
    fig, (ax0, ax2) = plt.subplots(1, 2)
    p = FigureWindow("test", custom_figure=(fig, ax0))
    new_ax = p._add_subplot_axes(p.ax, loc="top")
    assert get_axis_geometry(ax0) == (2, 1, 1)
    assert get_axis_geometry(new_ax) == (2, 1, 0)
    assert get_axis_geometry(ax2) == (1, 2, 1)
    assert len(p.fig.get_axes()) == 3


def test_add_subplot_multiple_rows_top():
    fig, (ax0, ax1) = plt.subplots(2, 1)
    p = FigureWindow("test", custom_figure=(fig, ax0))
    new_ax = p._add_subplot_axes(p.ax, loc="top")
    assert get_axis_geometry(ax0) == (2, 1, 1)
    assert get_axis_geometry(new_ax) == (2, 1, 0)
    assert get_axis_geometry(ax1) == (2, 1, 1)
    assert len(fig.get_axes()) == 3


def test_add_subplot_multiple_rows_top_middle():
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1)
    p = FigureWindow("test", custom_figure=(fig, ax1))
    new_ax = p._add_subplot_axes(p.ax, loc="top")
    assert get_axis_geometry(ax0) == (3, 1, 0)
    assert get_axis_geometry(ax1) == (2, 1, 1)
    assert get_axis_geometry(new_ax) == (2, 1, 0)
    assert get_axis_geometry(ax2) == (3, 1, 2)
    assert len(fig.get_axes()) == 4


def test_add_subplot_matrix_top():
    fig, axes = plt.subplots(3, 3)
    p = FigureWindow("test", custom_figure=(fig, axes[1, 1]))
    new_ax = p._add_subplot_axes(p.ax, loc="top")
    assert get_axis_geometry(axes[0, 2]) == (3, 3, 2)
    assert get_axis_geometry(axes[1, 2]) == (3, 3, 5)
    assert get_axis_geometry(p.ax) == (2, 1, 1)
    assert get_axis_geometry(new_ax) == (2, 1, 0)
    assert len(fig.get_axes()) == 10


# ========== TESTS FOR BOTTOM SUBPLOT ADDITION ==========
def test_add_subplot_singleplot_bottom():
    p = FigureWindow("test")
    new_ax = p._add_subplot_axes(p.ax, loc="bottom")
    assert get_axis_geometry(p.ax) == (2, 1, 0)
    assert get_axis_geometry(new_ax) == (2, 1, 1)
    assert len(p.fig.get_axes()) == 2


def test_add_subplot_multiple_columns_bottom_end():
    fig, (ax0, ax1) = plt.subplots(1, 2)
    p = FigureWindow("test", custom_figure=(fig, ax1))
    new_ax = p._add_subplot_axes(p.ax, loc="bottom")
    assert get_axis_geometry(ax0) == (1, 2, 0)
    assert get_axis_geometry(ax1) == (2, 1, 0)
    assert get_axis_geometry(new_ax) == (2, 1, 1)
    assert len(p.fig.get_axes()) == 3


def test_add_subplot_multiple_columns_bottom_front():
    fig, (ax0, ax2) = plt.subplots(1, 2)
    p = FigureWindow("test", custom_figure=(fig, ax0))
    new_ax = p._add_subplot_axes(p.ax, loc="bottom")
    assert get_axis_geometry(ax0) == (2, 1, 0)
    assert get_axis_geometry(new_ax) == (2, 1, 1)
    assert get_axis_geometry(ax2) == (1, 2, 1)
    assert len(p.fig.get_axes()) == 3


def test_add_subplot_multiple_rows_bottom():
    fig, (ax0, ax1) = plt.subplots(2, 1)
    p = FigureWindow("test", custom_figure=(fig, ax0))
    new_ax = p._add_subplot_axes(p.ax, loc="bottom")
    assert get_axis_geometry(ax0) == (2, 1, 0)
    assert get_axis_geometry(new_ax) == (2, 1, 1)
    assert get_axis_geometry(ax1) == (2, 1, 1)
    assert len(fig.get_axes()) == 3


def test_add_subplot_multiple_rows_bottom_middle():
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1)
    p = FigureWindow("test", custom_figure=(fig, ax1))
    new_ax = p._add_subplot_axes(p.ax, loc="bottom")
    assert get_axis_geometry(ax0) == (3, 1, 0)
    assert get_axis_geometry(ax1) == (2, 1, 0)
    assert get_axis_geometry(new_ax) == (2, 1, 1)
    assert get_axis_geometry(ax2) == (3, 1, 2)
    assert len(fig.get_axes()) == 4


def test_add_subplot_matrix_bottom():
    fig, axes = plt.subplots(3, 3)
    p = FigureWindow("test", custom_figure=(fig, axes[2, 1]))
    new_ax = p._add_subplot_axes(p.ax, loc="bottom")
    assert get_axis_geometry(axes[0, 2]) == (3, 3, 2)
    assert get_axis_geometry(axes[1, 2]) == (3, 3, 5)
    assert get_axis_geometry(p.ax) == (2, 1, 0)
    assert get_axis_geometry(new_ax) == (2, 1, 1)
    assert len(fig.get_axes()) == 10
