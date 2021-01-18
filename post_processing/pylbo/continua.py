import numpy as np
import matplotlib.collections as mpl_collections
import matplotlib.lines as mpl_lines
import matplotlib.patches as mpl_patches
from matplotlib import colors

from pylbo._version import _mpl_version
from pylbo.utilities.logger import pylboLogger

CONTINUA_NAMES = ["slow-", "slow+", "alfven-", "alfven+", "thermal", "doppler"]
CONTINUA_COLORS = ["red", "red", "cyan", "cyan", "green", "grey"]


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
        CONTINUA_NAMES[0]: doppler - np.sqrt(slow2),
        CONTINUA_NAMES[1]: doppler + np.sqrt(slow2),
        CONTINUA_NAMES[2]: doppler - np.sqrt(alfven2),
        CONTINUA_NAMES[3]: doppler + np.sqrt(alfven2),
        CONTINUA_NAMES[4]: thermal,
        CONTINUA_NAMES[5]: doppler,
    }
    return continua


class ContinuaHandler:
    """
    Handler to draw continua regions on the plots and make them interactive.

    Parameters
    ----------
    interactive : bool
        If `True`, makes the legend pickable and continuum plotting interactive.

    Attributes
    ----------
    legend : `matplotlib.legend.Legend`
        The matplotlib legend.
    legend_properties : dict
        Properties related to the legend.
    """

    def __init__(self, interactive):
        self.legend = None
        self.alpha_point = 0.8
        self.alpha_region = 0.2
        self.alpha_hidden = 0.05

        self.marker = "p"
        self.markersize = 64
        self.pickradius = 10
        self.linewidth = 2
        self.legend_properties = {}

        self.interactive = interactive
        self.continua_names = CONTINUA_NAMES
        self._continua_colors = CONTINUA_COLORS
        self._drawn_items = []
        self._legend_mapping = {}

    def on_legend_pick(self, event):
        """
        Determines what happens when the legend gets picked.

        Parameters
        ----------
        event : `~matplotlib.backend_bases.PickEvent`
            The matplotlib pick event.
        """
        artist = event.artist
        if artist not in self._legend_mapping:
            return
        drawn_item = self._legend_mapping.get(artist)
        visible = not drawn_item.get_visible()
        drawn_item.set_visible(visible)
        if visible:
            if isinstance(artist, (mpl_collections.PathCollection, mpl_lines.Line2D)):
                artist.set_alpha(self.alpha_point)
            else:
                artist.set_alpha(self.alpha_region)
        else:
            artist.set_alpha(self.alpha_hidden)
        artist.figure.canvas.draw()

    def make_legend_pickable(self):
        """
        Makes the legend pickable, only used if interactive.
        """
        legend_handles = self.legend.legendHandles
        handle_labels = [handle.get_label() for handle in legend_handles]
        # we need a mapping of the legend item to the actual item that was drawn
        for i, drawn_item in enumerate(self._drawn_items):
            # TODO: for some reason fill_between returns empty handles such that
            #       the code below errors out. The try-except clause is an attempt
            #       to fix it and _should_ work. This relies on the ordening of the
            #       legend items when they are drawn, but should be thorougly tested.
            try:
                idx = handle_labels.index(drawn_item.get_label())
                legend_item = self.legend.legendHandles[idx]
            except ValueError:
                idx = i
                legend_item = self.legend.legendHandles[idx]
                # fix empty label
                legend_item.set_label(drawn_item.get_label())
            if not drawn_item.get_label() == legend_item.get_label():
                raise ValueError(
                    f"something went wrong in mapping legend items to drawn items. \n"
                    f"Tried to map {legend_item} (label '{legend_item.get_label()}')"
                    f" to {drawn_item} (label '{drawn_item.get_label()}') \n"
                )
            # set_picker is deprecated for line2D from matplotlib 3.3 onwards
            if isinstance(legend_item, mpl_lines.Line2D) and _mpl_version >= "3.3":
                legend_item.set_picker(True)
                legend_item.pickradius = self.pickradius
            else:
                legend_item.set_picker(self.pickradius)
            legend_item.set_alpha(self.alpha_hidden)
            # we make the continuum regions invisible until clicked
            drawn_item.set_visible(False)
            self._legend_mapping[legend_item] = drawn_item

    def add(self, item):
        """
        Adds an item to the list of drawn items on the canvas.

        Parameters
        ----------
        item : object
            A single object, usually a return from the matplotlib plot or scatter
            methods.
        """
        if isinstance(item, (list, np.ndarray, tuple)):
            raise ValueError("object expected, not something list-like")
        self._drawn_items.append(item)

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
        continuum : np.ndarray
            Array with values.

        Returns
        -------
        True if all values are the same, false otherwise.
        """
        if all(np.diff(continuum) == 0):
            return True
        return False

    @staticmethod
    def check_if_all_zero(continuum):
        """
        Checks if a given continuum is pure zero.

        Parameters
        ----------
        continuum : np.ndarray
            Array with values.

        Returns
        -------
        True if all values are zero, false otherwise.
        """
        if all(continuum == 0):
            return True
        return False
