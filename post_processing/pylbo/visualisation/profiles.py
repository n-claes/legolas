import numpy as np
from pylbo.visualisation.figure_manager import FigureWindow
from pylbo.visualisation.legend_interface import LegendHandler
from pylbo.visualisation.continua import ContinuaHandler


class EquilibriumProfile(FigureWindow):
    """Subclass responsible for drawing the equilibrium profiles."""

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
                    self.ax2 = super()._add_subplot_axes(self.ax, "bottom")
                    break
        self.leg_handle = LegendHandler(interactive)
        self.draw()
        if interactive:
            self._enable_interactive_legend(self.leg_handle)

    def draw(self):
        """Draws the figure."""
        super().draw()
        self._add_equilibria()
        self.fig.tight_layout()

    def _add_equilibria(self):
        """
        Adds the equilibria to the figure. Also sets the legend handler items
        """
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
        # we add the plasma beta as well, only for MHD
        if not np.all(np.isclose(self.data.equilibria["B0"], 0)):
            plasma_beta = (
                2 * self.data.equilibria["rho0"] * self.data.equilibria["T0"]
            ) / self.data.equilibria["B0"] ** 2
            (item,) = self.ax.plot(
                self.data.grid_gauss,
                plasma_beta,
                label=r"plasma-$\beta$",
                alpha=self.leg_handle.alpha_point,
            )
            self.leg_handle.add(item)
            items.append(item)
        labels = [l_item.get_label() for l_item in items]
        self.leg_handle.legend = self.ax.legend(
            items,
            labels,
            bbox_to_anchor=(0.0, 1, 1, 0.102),
            loc="lower left",
            ncol=10,
            mode="expand",
        )
        # enable autoscaling when clicking
        self.leg_handle.autoscale = True


class ContinuumProfile(FigureWindow):
    """Subclass responsible for drawing the continuum profiles."""

    def __init__(self, data, figsize, interactive, **kwargs):
        super().__init__(figure_type="continua", figsize=figsize)
        self.data = data
        self.kwargs = kwargs
        self.handler = ContinuaHandler(interactive)
        self.draw()
        if interactive:
            self._enable_interactive_legend(self.handler)

    def draw(self):
        """Draws the continua."""
        super().draw()
        self._draw_continua()
        self.fig.tight_layout()

    def _draw_continua(self):
        """Adds the continua to the plot, also sets the legend handlers."""
        for color, name in zip(
            self.handler.continua_colors, self.handler.continua_names
        ):
            continuum = self.data.continua[name]
            if self.handler.check_if_all_zero(continuum):
                continue
            # non-adiabatic slow continua have real and imaginary parts
            if np.any(np.iscomplex(continuum)) and "slow" in name:
                (item,) = self.ax.plot(
                    self.data.grid_gauss,
                    self.data.continua[name].real,
                    color=color,
                    label="".join([name, r"$_{Re}$"]),
                )
                self.handler.add(item)
                (item,) = self.ax.plot(
                    self.data.grid_gauss,
                    self.data.continua[name].imag,
                    linestyle="dashed",
                    color=color,
                    label="".join([name, r"$_{Im}$"]),
                )
                self.handler.add(item)
            else:
                cont = self.data.continua[name]
                if name == "thermal":
                    cont = cont.imag
                (item,) = self.ax.plot(
                    self.data.grid_gauss, cont, color=color, label=name
                )
                self.handler.add(item)
        self.handler.legend = self.ax.legend()
        self.handler.autoscale = True
        self.ax.axhline(y=0, color="grey", linestyle="dotted", alpha=0.8)
        self.ax.axvline(
            x=self.data.x_start, color="grey", linestyle="dotted", alpha=0.8
        )
        self.ax.set_xlabel("Grid coordinate")
        self.ax.set_ylabel(r"$\omega$")


class EquilibriumBalance(FigureWindow):
    """Subclass responsible for plotting the equilibrium balance equations."""

    def __init__(self, data, figsize, **kwargs):
        super().__init__(figure_type="equilibrium-balance", figsize=figsize)
        self.data = data
        self.kwargs = kwargs
        self.ax2 = super()._add_subplot_axes(self.ax, "bottom")
        self.draw()

    def draw(self):
        """Draws the equilibrium balance equations."""
        rho = self.data.equilibria["rho0"]
        drho = self.data.equilibria["drho0"]
        temp = self.data.equilibria["T0"]
        dtemp = self.data.equilibria["dT0"]
        b02 = self.data.equilibria["B02"]
        db02 = self.data.equilibria["dB02"]
        b03 = self.data.equilibria["B03"]
        db03 = self.data.equilibria["dB03"]
        g = self.data.equilibria["grav"]
        v02 = self.data.equilibria["v02"]
        kappa_perp = self.data.equilibria["kappa_perp"]
        # L0 is only non-zero when custom heating is added
        # (and is not saved to the datfile for now)
        heat_loss = np.zeros_like(self.data.grid_gauss)
        r_scale = self.data.scale_factor
        dr_scale = self.data.d_scale_factor
        equil_force = (
            drho * temp
            + rho * dtemp
            + b02 * db02
            + b03 * db03
            + rho * g
            - (dr_scale / r_scale) * (rho * v02 ** 2 - b02 ** 2)
        )
        equil_force[np.where(abs(equil_force) <= 1e-16)] = 0
        self.ax.plot(self.data.grid_gauss, equil_force, **self.kwargs)
        self.ax.axhline(y=0, color="grey", linestyle="dotted")
        if any(abs(equil_force) > 1e-14):
            self.ax.set_yscale("symlog")
        self.ax.set_title("Force balance")

        # ddT0 is not saved, so we do it numerically (it's a check anyway)
        dtemp_fact = np.gradient(kappa_perp * dtemp, self.data.grid_gauss, edge_order=2)
        equil_nadiab = (
            dr_scale * kappa_perp * dtemp / r_scale + dtemp_fact - rho * heat_loss
        )
        equil_nadiab[np.where(abs(equil_nadiab) <= 1e-16)] = 0
        self.ax2.plot(self.data.grid_gauss, equil_nadiab, **self.kwargs)
        self.ax2.axhline(y=0, color="grey", linestyle="dotted")
        if any(abs(equil_nadiab) > 1e-14):
            self.ax2.set_yscale("symlog")
        self.ax2.set_title("Nonadiabatic balance")
