from pylbo.visualisation.figure_manager import FigureWindow
from pylbo.visualisation.legend_interface import LegendHandler


class EquilibriumProfile(FigureWindow):
    def __init__(self, data, figsize, interactive, **kwargs):
        super().__init__(figure_type="equilibrium-fields", figsize=figsize)
        self.data = data
        self.kwargs = kwargs
        self.ax2 = None
        # check if we need an additional axis
        for name in self.data.eq_names:
            # check derivative names
            if name.startswith("d"):
                values = self.data.equilibria.get(name)
                if any(values != 0):
                    self.ax2 = self._add_extra_axis()
                    break
        self.leg_handle = LegendHandler(interactive)
        self.draw()
        if interactive:
            self._enable_interactive_legend()

    def _add_extra_axis(self):
        self.ax.change_geometry(2, 1, 1)
        ax2 = self.fig.add_subplot(212)
        self.fig.tight_layout()
        return ax2

    def draw(self):
        super().draw()
        self._add_equilibria()

    def _add_equilibria(self):
        items = []
        for name in self.data.eq_names:
            if name.startswith("d"):
                axis = self.ax2
            else:
                axis = self.ax
            values = self.data.equilibria.get(name)
            if all(values == 0):
                continue
            (item,) = axis.plot(
                self.data.grid_gauss,
                values,
                label=name,
                alpha=self.leg_handle.alpha_point,
            )
            axis.axhline(y=0, color="grey", linestyle="dotted", alpha=0.6)
            axis.axvline(
                x=self.data.x_start, color="grey", linestyle="dotted", alpha=0.6
            )
            self.leg_handle.add(item)
            items.append(item)
        labels = [l_item.get_label() for l_item in items]
        self.leg_handle.legend = self.ax.legend(
            items,
            labels,
            bbox_to_anchor=(0., 1, 1, .102),
            loc="lower left",
            ncol=10,
            mode="expand",
        )
        # enable autoscaling when clicking
        self.leg_handle.autoscale = True
        self.fig.tight_layout()

    def _enable_interactive_legend(self):
        self.leg_handle.make_legend_pickable()
        callback_kind = "pick_event"
        callback_method = self.leg_handle.on_legend_pick
        callback_id = self.fig.canvas.mpl_connect(callback_kind, callback_method)
        self._mpl_callbacks.append(
            {
                "cid": callback_id,
                "kind": callback_kind,
                "method": callback_method,
            }
        )
