import numpy as np
from pylbo.utilities.logger import pylboLogger

CONTINUA_KEYS = ["slow-", "slow+", "alfven-", "alfven+", "thermal", "doppler"]
CONTINUA_COLOURS = ["red", "red", "cyan", "cyan", "green", "grey"]


def calculate_continua(ds):
    """
    Calculates the different continua for a given dataset.

    Parameters
    ----------
    ds : `~pylbo.data_containers.LegolasDataSet` instance
        The Legolas dataset.

    Returns
    -------
    continua : dict
        Dictonary containing the various continua as numpy arrays.
    """
    rho = ds.equilibria["rho0"]
    B02 = ds.equilibria["B02"]
    B03 = ds.equilibria["B03"]
    B0 = np.sqrt(B02 ** 2 + B03 ** 2)
    v02 = ds.equilibria["v02"]
    v03 = ds.equilibria["v03"]
    T = ds.equilibria["T0"]
    p = rho * T

    k2 = ds.parameters["k2"]
    k3 = ds.parameters["k3"]
    gamma = ds.gamma

    # Alfven and slow continuum (squared)
    alfven2 = (1 / rho) * ((k2 * B02 / ds.scale_factor) + k3 * B03) ** 2
    slow2 = (gamma * p / (gamma * p + B0 ** 2)) * alfven2
    # doppler shift equals dot product of k and v
    doppler = k2 * v02 / ds.scale_factor + k3 * v03

    # Thermal continuum
    dLdT = ds.equilibria["dLdT"]
    dLdrho = ds.equilibria["dLdrho"]
    kappa_para = ds.equilibria["kappa_para"]
    kappa_perp = ds.equilibria["kappa_perp"]
    # if there is no conduction/cooling, there is no thermal continuum.
    if (
        all(dLdT == 0)
        and all(dLdrho == 0)
        and all(kappa_para == 0)
        and all(kappa_perp == 0)
    ):
        thermal = np.zeros_like(ds.grid_gauss)
    # if temperature is zero (no pressure), set to zero and return
    elif (T == 0).all():
        thermal = np.zeros_like(ds.grid_gauss)
    else:
        # wave vector parallel to magnetic field,
        # uses vector projection and scale factor
        kpara = (k2 * B02 / ds.scale_factor + k3 * B03) / B0
        cs2 = gamma * p / rho  # sound speed
        vA2 = B0 ** 2 / rho  # Alfvén speed
        ci2 = p / rho  # isothermal sound speed
        sigma_A2 = kpara * vA2  # Alfvén frequency
        sigma_c2 = cs2 * sigma_A2 / (vA2 + cs2)  # cusp frequency
        sigma_i2 = ci2 * sigma_A2 / (vA2 + ci2)  # isothermal cusp frequency

        # thermal and slow continuum are coupled in a third degree polynomial in omega,
        # coeffi means the coefficient corresponding to the term omega^i
        coeff3 = rho * (cs2 + vA2) * 1j / (gamma - 1)
        coeff2 = (
            -(kappa_para * kpara ** 2 + rho * dLdT) * (ci2 + vA2) + rho ** 2 * dLdrho
        )
        coeff1 = -rho * (cs2 + vA2) * sigma_c2 * 1j / (gamma - 1)
        coeff0 = (kappa_para * kpara ** 2 + rho * dLdT) * (
            ci2 + vA2
        ) * sigma_i2 - rho ** 2 * dLdrho * sigma_A2
        # we have to solve this equation "gauss_gridpts" times.
        # the thermal continuum corresponds to the (only) purely imaginary solution,
        # slow continuum are other two (real) solutions
        thermal = []
        for idx in range(len(ds.grid_gauss)):
            solutions = np.roots([coeff3[idx], coeff2[idx], coeff1[idx], coeff0[idx]])
            imag_sol = solutions.imag[abs(solutions.real) < 1e-14]
            if (imag_sol == 0).all():
                # if all solutions are zero (no thermal continuum), then append zero
                w = 0
            else:
                # else there should be ONE value that is nonzero,
                # this is the thermal continuum value.
                # Filter out non-zero solution and try to unpack.
                # If there is more than one non-zero
                # value unpacking fails, this is a sanity check so we raise an error.
                try:
                    (w,) = imag_sol[imag_sol != 0]  # this is an array, so unpack with ,
                except ValueError:
                    pylboLogger.warning(
                        f"Something went wrong, more than one solution for the thermal "
                        f"continuum was found (and there should be only one). "
                        f"Setting value to zero. "
                        f"Solutions found: {imag_sol[imag_sol != 0]}."
                    )
                    w = 0
            thermal.append(w)
        thermal = np.asarray(thermal)

        # additional sanity check: if the slow continuum vanishes
        # there is an analytical solution for the
        # thermal continuum. This one should be equal to the solution
        # obtained through solving the polynomial equation.
        if (slow2 == 0).all() and (T != 0).all():
            analytic_sol = (
                1j
                * (gamma - 1)
                * (rho * dLdrho / T - dLdT * (T + vA2))
                / (vA2 + gamma * T)
            )
            diff = abs(np.sort(thermal) - np.sort(np.imag(analytic_sol)))
            # if any value in the difference is larger than 1e-12, print a warning
            if any(diff > 1e-12):
                pylboLogger.warning(
                    "slow continuum vanishes, but difference in analytical/numerical "
                    "thermal continuum!"
                )

    # get doppler-shifted continua and return
    continua = {
        CONTINUA_KEYS[0]: doppler - np.sqrt(slow2),
        CONTINUA_KEYS[1]: doppler + np.sqrt(slow2),
        CONTINUA_KEYS[2]: doppler - np.sqrt(alfven2),
        CONTINUA_KEYS[3]: doppler + np.sqrt(alfven2),
        CONTINUA_KEYS[4]: thermal,
        CONTINUA_KEYS[5]: doppler,
    }
    return continua


class ContinuaHandler:
    def __init__(self, figure_obj):
        self.figure_obj = figure_obj
        self.legend_items = []
        self.legend = None
        self.key_colours = CONTINUA_COLOURS

    def on_legend_pick(self, event):
        # TODO: add interactivity with legend items
        pass

    def draw_continua(self, **kwargs):
        for key, colour in zip(CONTINUA_KEYS, self.continua_colours):
            continuum = self.figure_obj.dataset.continua[key]
            if self._check_if_all_zero(continuum=continuum):
                continue
            min_value = np.min(continuum)
            max_value = np.max(continuum)
            if self._check_if_collapsed(continuum=continuum):
                item = self.figure_obj.ax.plot(
                    min_value,
                    0,
                    marker=kwargs.pop("marker", "p"),
                    markersize=kwargs.pop("markersize", 8),
                    color=colour,
                    linestyle="none",
                    alpha=kwargs.pop("alpha", 0.5),
                    label=key,
                )
                self.legend_items.append(item)
            else:
                raise NotImplementedError
            self.legend = self.figure_obj.ax.legend(loc=kwargs.pop("legend_loc", "best"))

    @property
    def continua_colours(self):
        return self.key_colours

    @continua_colours.setter
    def continua_colours(self, values):
        if not isinstance(values, (list, np.ndarray)):
            raise ValueError(
                f"continua_colours should be an array/list but got {type(values)}"
            )
        self.continua_colours = values

    @staticmethod
    def _check_if_collapsed(continuum):
        if all(np.diff(continuum) == 0):
            return True
        return False

    @staticmethod
    def _check_if_all_zero(continuum):
        if all(continuum == 0):
            return True
        return False
