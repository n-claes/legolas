from __future__ import annotations

import difflib
from typing import Union

import numpy as np
from pylbo.data_containers import LegolasDataSet
from pylbo.utilities.logger import pylboLogger
from pylbo.visualisation.utils import ef_name_to_latex, validate_ef_name


class ModeVisualisationData:
    """
    Class that contains the data used for eigenmode visualisations.

    Parameters
    ----------
    ds : ~pylbo.data_containers.LegolasDataSet
        The dataset containing the eigenfunctions and modes to visualise.
    omega : list[complex]
        The (approximate) eigenvalue(s) of the mode(s) to visualise.
    ef_name : str
        The name of the eigenfunction to visualise.
    use_real_part : bool
        Whether to use the real part of the eigenmode solution.
    complex_factor : complex
        A complex factor to multiply the eigenmode solution with.
    add_background : bool
        Whether to add the equilibrium background to the eigenmode solution.

    Attributes
    ----------
    ds : ~pylbo.data_containers.LegolasDataSet
        The dataset containing the eigenfunctions and modes to visualise.
    omega : list[complex]
        The (approximate) eigenvalue(s) of the mode(s) to visualise.
    eigenfunction : list[np.ndarray]
        The eigenfunction of the mode(s) to visualise.
    use_real_part : bool
        Whether to use the real part of the eigenmode solution.
    complex_factor : complex
        The complex factor to multiply the eigenmode solution with.
    add_background : bool
        Whether to add the equilibrium background to the eigenmode solution.
    """

    def __init__(
        self,
        ds: LegolasDataSet,
        omega: list[complex],
        ef_name: str = None,
        use_real_part: bool = True,
        complex_factor: complex = None,
        add_background: bool = False,
    ) -> None:
        self.ds = ds
        self.use_real_part = use_real_part
        self.complex_factor = self._validate_complex_factor(complex_factor)
        self.add_background = add_background
        self._print_bg_info = True

        self._ef_name = None if ef_name is None else validate_ef_name(ds, ef_name)
        self._ef_name_latex = None if ef_name is None else self.get_ef_name_latex()
        self._all_efs = self._get_all_efs(ds, omega)
        self.omega = [all_efs.get("eigenvalue") for all_efs in self._all_efs]
        self.eigenfunction = [all_efs.get(self._ef_name) for all_efs in self._all_efs]

    @property
    def k2(self) -> float:
        """The k2 wave number of the eigenmode solution."""
        return self.ds.parameters["k2"]

    @property
    def k3(self) -> float:
        """The k3 wave number of the eigenmode solution."""
        return self.ds.parameters["k3"]

    @property
    def part_name(self) -> str:
        """
        Returns the name of the part of the eigenmode solution to use, i.e.
        'real' or 'imag'.
        """
        return "real" if self.use_real_part else "imag"

    def _get_all_efs(self, ds, omega):
        """
        Returns an array of dicts with all eigenfunctions for every eigenvalue.
        The dictionaries will be updated with the derived eigenfunctions if they
        are available in the dataset.

        Parameters
        ----------
        ds : ~pylbo.data_containers.LegolasDataSet
            The dataset containing the eigenfunctions.
        omega : list[complex]
            The (approximate) eigenvalue(s) of the mode(s) to retrieve the
            eigenfunctions from.

        Returns
        -------
        np.ndarray
            An array of dicts with all eigenfunctions for every eigenvalue.
        """
        arr1 = ds.get_eigenfunctions(omega)
        if not ds.derived_efs_written:
            return arr1

        arr2 = ds.get_derived_eigenfunctions(omega)
        arr = np.empty(len(arr1) + len(arr2), dtype=dict)
        for i, (dict1, dict2) in enumerate(zip(arr1, arr2)):
            ev1 = dict1.get("eigenvalue")
            ev2 = dict2.get("eigenvalue")
            if not np.isclose(ev1, ev2, atol=1e-12):
                pylboLogger.warning(
                    f"The eigenvalue of the eigenfunction {ev1:.6e} and the derived "
                    f"eigenfunction {ev2:.6e} do not match. Using eigenfunctions only."
                )
                return arr1
            arr[i] = {**dict1, **dict2}
        return arr

    def get_ef_name_latex(self) -> str:
        """Returns the latex representation of the eigenfunction name."""
        return ef_name_to_latex(
            self._ef_name, geometry=self.ds.geometry, real_part=self.use_real_part
        )

    def _validate_complex_factor(self, complex_factor: complex) -> complex:
        """
        Validates the complex factor.

        Parameters
        ----------
        complex_factor : complex
            The complex factor to validate.

        Returns
        -------
        complex
            The complex factor if it is valid, otherwise 1.
        """
        return complex_factor if complex_factor is not None else 1

    def get_mode_solution(
        self,
        ef: np.ndarray,
        omega: complex,
        u2: Union[float, np.ndarray],
        u3: Union[float, np.ndarray],
        t: Union[float, np.ndarray],
    ) -> np.ndarray:
        """
        Calculates the full eigenmode solution for given coordinates and time.
        If a complex factor was given, the eigenmode solution is multiplied with the
        complex factor. If :attr:`use_real_part` is True the real part of the eigenmode
        solution is returned, otherwise the complex part. If :attr:`add_background` is
        True, the background is added to the eigenmode solution.

        Parameters
        ----------
        ef : np.ndarray
            The eigenfunction to use.
        omega : complex
            The eigenvalue to use.
        u2 : Union[float, np.ndarray]
            The y coordinate(s) of the eigenmode solution.
        u3 : Union[float, np.ndarray]
            The z coordinate(s) of the eigenmode solution.
        t : Union[float, np.ndarray]
            The time(s) of the eigenmode solution.

        Returns
        -------
        np.ndarray
            The real or imaginary part of the eigenmode solution for the given
            set of coordinate(s) and time(s).
        """
        solution = (
            self.complex_factor
            * ef
            * np.exp(1j * self.k2 * u2 + 1j * self.k3 * u3 - 1j * omega * t)
        )
        return getattr(solution, self.part_name)

    def get_background(self, shape: tuple[int, ...], name=None) -> np.ndarray:
        """
        Returns the background of the eigenmode solution.

        Parameters
        ----------
        shape : tuple[int, ...]
            The shape of the eigenmode solution.
        name : str
            The name of the background to use. If None, the background name
            will be inferred from the eigenfunction name.

        Returns
        -------
        np.ndarray
            The background of the eigenmode solution, sampled on the eigenfunction
            grid and broadcasted to the same shape as the eigenmode solution.
        """
        if name is None:
            name = self._get_background_name()
        bg = self.ds.equilibria[name]
        bg_sampled = self._sample_background_on_ef_grid(bg)
        if self._print_bg_info:
            pylboLogger.info(f"background {name} broadcasted to shape {shape}")
        return np.broadcast_to(bg_sampled, shape=reversed(shape)).transpose()

    def _sample_background_on_ef_grid(self, bg: np.ndarray) -> np.ndarray:
        """
        Samples the background array on the eigenfunction grid.

        Parameters
        ----------
        bg : np.ndarray
            The background array with Gaussian grid spacing

        Returns
        -------
        np.ndarray
            The background array with eigenfunction grid spacing
        """
        if self._print_bg_info:
            pylboLogger.info(
                f"sampling background [{len(bg)}] on eigenfunction grid "
                f"[{len(self.ds.ef_grid)}]"
            )
        return np.interp(self.ds.ef_grid, self.ds.grid_gauss, bg)

    def _get_background_name(self) -> str:
        """
        Returns the name of the background.

        Returns
        -------
        str
            The closest match between the eigenfunction name and the equilibrium
            name.

        Raises
        ------
        ValueError
            If the eigenfunction name is a magnetic vector potential component.
        """
        if self._ef_name in ("a1", "a2", "a3"):
            raise ValueError(
                "Unable to add a background to the magnetic vector potential."
            )
        (name,) = difflib.get_close_matches(self._ef_name, self.ds.eq_names, 1)
        if self._print_bg_info:
            pylboLogger.info(
                f"adding background for '{self._ef_name}', closest match is '{name}'"
            )
        return name
