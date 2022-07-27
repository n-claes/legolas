from typing import Union

import numpy as np
from pylbo.data_containers import LegolasDataSet
from pylbo.visualisation.utils import ef_name_to_latex


class ModeVisualisationData:
    """
    Class that contains the data used for eigenmode visualisations.

    Parameters
    ----------
    ds : ~pylbo.data_containers.LegolasDataSet
        The dataset containing the eigenfunctions and modes to visualise.
    omega : complex
        The (approximate) eigenvalue of the mode to visualise.
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
    omega : complex
        The (approximate) eigenvalue of the mode to visualise.
    eigenfunction : np.ndarray
        The eigenfunction to visualise.
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
        omega: complex,
        ef_name: str,
        use_real_part: bool = True,
        complex_factor: complex = None,
        add_background: bool = False,
    ) -> None:
        self.ds = ds
        self.omega, self.eigenfunction = self._retrieve_eigenfunction(omega, ef_name)
        self.use_real_part = use_real_part
        self.complex_factor = self._validate_complex_factor(complex_factor)
        self.add_background = add_background

        self._ef_name = ef_name
        self._ef_name_latex = self.get_ef_name_latex()

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

    def get_ef_name_latex(self) -> str:
        """
        Returns the latex representation of the eigenfunction name.
        """
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

    def _retrieve_eigenfunction(
        self, omega: complex, ef_name: str
    ) -> tuple[complex, np.ndarray]:
        """
        Retrieve the eigenfunction from the dataset.

        Parameters
        ----------
        omega : complex
            The (approximate) eigenvalue of the mode to visualise.
        ef_name : str
            The name of the eigenfunction to visualise.

        Returns
        -------
        complex
            The eigenvalue of the mode to visualise.
        np.ndarray
            The eigenfunction to visualise.
        """
        (efs,) = self.ds.get_eigenfunctions(omega)
        return efs.get("eigenvalue"), efs.get(ef_name)

    def get_mode_solution(
        self,
        ef: np.ndarray,
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
            set of coorinate(s) and time(s).
        """
        solution = (
            self.complex_factor
            * ef
            * np.exp(1j * self.k2 * u2 + 1j * self.k3 * u3 - 1j * self.omega * t)
        )
        if self.add_background:
            raise NotImplementedError()
        return getattr(solution, self.part_name)
