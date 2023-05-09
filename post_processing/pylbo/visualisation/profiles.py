import numpy as np
from pylbo.utilities.logger import pylboLogger
from pylbo.visualisation.continua import ContinuaHandler
from pylbo.visualisation.figure_window import FigureWindow, InteractiveFigureWindow
from pylbo.visualisation.legend_handler import LegendHandler
from pylbo.visualisation.utils import background_name_to_latex


class EquilibriumProfile(InteractiveFigureWindow):
    """Subclass responsible for drawing the equilibrium profiles."""

    def __init__(self, data, figsize, interactive, **kwargs):
        fig, ax = super().create_default_figure(
            figlabel="equilibrium-fields", figsize=figsize
        )
        super().__init__(fig)
        self.data = data
        self.kwargs = kwargs
        self.ax = ax
        self.ax2 = None
        # check if we need an additional axis
        for name in self.data.eq_names:
            # check derivative names
            if name.startswith("d"):
                values = self.data.equilibria.get(name)
                if any(values != 0):
                    self.ax2 = super().add_subplot_axes(self.ax, "bottom")
                    break
        self.leg_handle = LegendHandler(interactive)
        self.draw()
        if interactive:
            super().make_legend_interactive(self.leg_handle)
        self.fig.tight_layout()

    def draw(self):
        """Adds the equilibria to the figure. Also sets the legend handler items"""
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
                label=background_name_to_latex(name),
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


class ContinuumProfile(InteractiveFigureWindow):
    """Subclass responsible for drawing the continuum profiles."""

    def __init__(self, data, figsize, interactive, **kwargs):
        fig, ax = super().create_default_figure(figlabel="continua", figsize=figsize)
        super().__init__(fig)
        self.ax = ax
        self.data = data
        self.kwargs = kwargs
        self.handler = ContinuaHandler(interactive)
        self.draw()
        if interactive:
            self.make_legend_interactive(self.handler)
        self.fig.tight_layout()

    def draw(self):
        """Adds the continua to the plot, also sets the legend handlers."""
        for color, name in zip(
            self.handler.continua_colors, self.handler.continua_names
        ):
            continuum = self.data.continua[name]
            if np.allclose(abs(continuum), 0, atol=1e-12):
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
        fig, ax = super().create_default_figure(
            figlabel="equilibrium-balance", figsize=figsize
        )
        super().__init__(fig)
        self.data = data
        self.kwargs = kwargs
        self.ax = ax
        self.ax2 = super().add_subplot_axes(self.ax, "bottom")
        self.draw()

    def draw(self):
        """Draws the equilibrium balance equations."""
        force_balance = self._get_force_balance()
        self.ax.plot(self.data.grid_gauss, force_balance, **self.kwargs)
        self.ax.axhline(y=0, color="grey", linestyle="dotted")
        if any(abs(force_balance) > 1e-14):
            self.ax.set_yscale("symlog")
        self.ax.set_title("Force balance")

        energy_balance = self._get_energy_balance()
        self.ax2.plot(self.data.grid_gauss, energy_balance, **self.kwargs)
        self.ax2.axhline(y=0, color="grey", linestyle="dotted")
        if any(abs(energy_balance) > 1e-14):
            self.ax2.set_yscale("symlog")
        self.ax2.set_title("Energy (thermal) balance")

    def _get_force_balance(self):
        background = self.data.equilibria
        rho0 = background["rho0"]
        drho0 = background["drho0"]
        T0 = background["T0"]
        dT0 = background["dT0"]
        B02 = background["B02"]
        dB02 = background["dB02"]
        B03 = background["B03"]
        dB03 = background["dB03"]
        g0 = background["gravity"]
        v01 = background["v01"]
        dv01 = background["dv01"]
        v02 = background["v02"]
        eps = self.data.scale_factor
        deps = self.data.d_scale_factor
        return (
            drho0 * T0
            + rho0 * dT0
            + B02 * dB02
            + B03 * dB03
            + rho0 * g0
            - (deps / eps) * (rho0 * v02**2 - B02**2)
            + rho0 * v01 * dv01
        )

    def _get_energy_balance(self):
        background = self.data.equilibria
        rho0 = background["rho0"]
        T0 = background["T0"]
        dT0 = background["dT0"]
        ddT0 = background["ddT0"]
        v01 = background["v01"]
        dv01 = background["dv01"]
        B01 = background["B01"]
        B0 = background["B0"]
        L0 = background["L0"]
        tc_para = background["kappa_para"]
        tc_perp = background["kappa_perp"]

        dtc_perp_dr = background.get("dkappa_perp_dr", None)
        if dtc_perp_dr is None:
            pylboLogger.warning(
                "energy balance: dkappa_perp_dr not found in datfile. "
                "Deriving numerically using np.gradient."
            )
            dtc_perp_dr = np.gradient(tc_perp, self.data.grid_gauss, edge_order=2)

        is_mhd = "mhd" in self.data.header.get("physics_type", None)
        Kp = (tc_para - tc_perp) / B0**2 if is_mhd else np.zeros_like(tc_para)
        dKp = np.gradient(Kp, self.data.grid_gauss, edge_order=2)
        eps = self.data.scale_factor
        deps = self.data.d_scale_factor
        return (
            T0 * rho0 * (deps * v01 + eps * dv01) / eps
            + rho0 * L0
            - B01**2 * (Kp * dT0 + dKp * T0)
            - (1 / eps)
            * (deps * tc_perp * dT0 + eps * dtc_perp_dr * dT0 + eps * tc_perp * ddT0)
            + (1 / (self.data.gamma - 1)) * dT0 * rho0 * v01
        )
