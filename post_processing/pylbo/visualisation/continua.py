import numpy as np
from pylbo.utilities.logger import pylboLogger
from pylbo.visualisation.legend_interface import LegendHandler

CONTINUA_NAMES = ["slow-", "slow+", "alfven-", "alfven+", "thermal", "doppler"]
CONTINUA_COLORS = ["red", "red", "cyan", "cyan", "green", "grey"]


def calculate_continua(ds):
    """
    Calculates the different continua for a given dataset.

    Parameters
    ----------
    ds : ~pylbo.data_containers.LegolasDataSet
        The Legolas dataset.

    Returns
    -------
    continua : dict
        Dictonary containing the various continua as numpy arrays.
    """
    rho = ds.equilibria["rho0"]
    B02 = ds.equilibria["B02"]
    B03 = ds.equilibria["B03"]
    B0 = ds.equilibria["B0"]
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
    slow_min = -np.sqrt(slow2)
    slow_plus = np.sqrt(slow2)
    slow_is_zero = np.all(np.isclose(slow2, 0))

    # Thermal continuum calculation
    # wave vector parallel to magnetic field, uses vector projection and scale factor
    if not all(np.isclose(B0, 0, atol=1e-10)):
        # take care of hydro cases, this will throw a division by 0 warning
        kpara = (k2 * B02 / ds.scale_factor + k3 * B03) / B0
    else:
        kpara = 0
    cs2 = gamma * p / rho  # sound speed
    vA2 = B0 ** 2 / rho  # Alfvén speed
    ci2 = p / rho  # isothermal sound speed
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
    # if slow continuum vanishes, thermal continuum is analytical
    elif slow_is_zero:
        thermal = 1j * (gamma - 1) * (rho * dLdrho - dLdT * (ci2 + vA2)) / (cs2 + vA2)
    else:
        sigma_A2 = kpara ** 2 * vA2  # Alfvén frequency
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
        # modified slow continuum are other two solutions
        thermal = np.empty(shape=len(ds.grid_gauss), dtype=complex)
        slow_min = np.empty_like(thermal)
        slow_plus = np.empty_like(thermal)
        warning_logged = False
        for idx, _ in enumerate(ds.grid_gauss):
            solutions = np.roots([coeff3[idx], coeff2[idx], coeff1[idx], coeff0[idx]])
            # create mask for purely imaginary solutions
            mask = np.isclose(abs(solutions.real), 0)
            # this should only have 1 True element, the thermal solution
            if np.count_nonzero(mask) == 1:
                s_neg, s_pos = np.sort_complex(solutions[np.invert(mask)])
            else:
                if not warning_logged:
                    pylboLogger.warning(
                        f"encountered at least 1 point where the slow continuum has a "
                        f"real value close to zero. \nContinuum solutions: {solutions}"
                    )
                    warning_logged = True
                mask = np.array([False] * 3, dtype=bool)
                # check if this already happens on the first index
                if idx == 0:
                    # if so, assume largest mode is the thermal one and complex sort
                    # slow+ and slow- solutions
                    mask[solutions.imag.argmax()] = True
                    pylboLogger.warning(
                        f"this happened already at the first gridpoint, so we assumed "
                        f"{solutions[mask]} to be the thermal mode"
                    )
                    s_neg, s_pos = np.sort_complex(solutions[np.invert(mask)])
                else:
                    # if not, get solution closest to previous thermal mode
                    mask[np.abs(solutions - thermal[idx - 1]).argmin()] = True
                    # retrieve closest slow- solution
                    slow_solutions = solutions[np.invert(mask)]
                    idx_neg = np.abs(slow_solutions - slow_min[idx - 1]).argmin()
                    s_neg = slow_solutions[idx_neg]
                    # remaining value is slow+ solution, remove slow- from array
                    (s_pos,) = np.delete(slow_solutions, idx_neg)
            # extract continuum solutions
            (thermal[idx],) = solutions[mask]
            slow_min[idx] = s_neg
            slow_plus[idx] = s_pos

    # get doppler-shifted continua and return
    continua = {
        CONTINUA_NAMES[0]: doppler + slow_min,  # minus is accounted for in slow_min
        CONTINUA_NAMES[1]: doppler + slow_plus,
        CONTINUA_NAMES[2]: doppler - np.sqrt(alfven2),
        CONTINUA_NAMES[3]: doppler + np.sqrt(alfven2),
        CONTINUA_NAMES[4]: thermal,
        CONTINUA_NAMES[5]: doppler,
    }
    return continua


class ContinuaHandler(LegendHandler):
    """
    Handler to draw continua regions on the plots and make them interactive.

    Parameters
    ----------
    interactive : bool
        If `True`, makes the legend pickable and continuum plotting interactive.

    Attributes
    ----------
    continua_names : list
        The list of continua names
    """

    def __init__(self, interactive):
        super().__init__(interactive)
        self.continua_names = CONTINUA_NAMES
        self._continua_colors = CONTINUA_COLORS

    @property
    def continua_colors(self):
        """
        Returns the list of continua colors.

        Returns
        -------
        The continua colors as a list.
        """
        return self._continua_colors

    @continua_colors.setter
    def continua_colors(self, colors):
        """
        Setter for the continua colors attribute.

        Parameters
        ----------
        colors : list, numpy.ndarray
            The colors to use when plotting the continua as a list of strings.

        Raises
        ------
        ValueError
            If a wrong argument is passed or if it is of improper length.
        """
        if colors is None:
            return
        if not isinstance(colors, (list, np.ndarray)):
            raise ValueError(
                f"continua_colors should be an array/list but got {type(colors)}"
            )
        if not len(colors) == len(CONTINUA_COLORS):
            raise ValueError(
                f"continua_colors should be of length {len(CONTINUA_COLORS)}"
            )
        self._continua_colors = colors

    @staticmethod
    def check_if_collapsed(continuum):
        """
        Checks if a given continuum is "collapsed" to a single point.

        Parameters
        ----------
        continuum : numpy.ndarray
            Array with values.

        Returns
        -------
        True if all values are the same, false otherwise.
        """
        if np.all(np.isclose(np.diff(continuum), 0)):
            return True
        return False

    @staticmethod
    def check_if_all_zero(continuum):
        """
        Checks if a given continuum is pure zero.

        Parameters
        ----------
        continuum : numpy.ndarray
            Array with values.

        Returns
        -------
        True if all values are zero, false otherwise.
        """
        if all(continuum == 0):
            return True
        return False
