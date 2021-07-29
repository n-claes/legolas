import numpy as np
from pylbo.utilities.toolbox import transform_to_numpy
from pylbo.data_containers import LegolasDataSet
from pylbo.exceptions import EigenfunctionsNotPresent
from pylbo.visualisation.eigenfunction_interface import EigenfunctionInterface


class EigenfunctionHandler(EigenfunctionInterface):
    """
    Main handler for eigenfunctions.
    """

    def __init__(self, data, ef_ax):
        super().__init__(data, ef_ax)
        self._function_names = self.data.ef_names

    def _check_data_is_present(self):
        if not any(transform_to_numpy(self.data.efs_written)):
            raise EigenfunctionsNotPresent(
                "None of the given datfiles has eigenfunctions written to it."
            )

    def update_plot(self):
        self.axis.clear()
        if not self._selected_idxs:
            return
        ef_name = self._function_names[self._selected_name_idx]
        for ds, idxs_dict in self._selected_idxs.items():
            idxs = np.array([int(idx) for idx in idxs_dict.keys()])
            ef_container = ds.get_eigenfunctions(ev_idxs=idxs)
            for ev_idx, efs in zip(idxs, ef_container):
                ef = efs.get(ef_name)
                if self._use_real_part:
                    ef = ef.real
                else:
                    ef = ef.imag
                # check for retransform
                if (
                    self._retransform
                    and ef_name in ("rho", "v1", "v3", "T", "a2")
                    and self.data.geometry == "cylindrical"
                ):
                    ef = ef * ds.ef_grid
                label = super()._get_label(ds, ev_idx, efs.get("eigenvalue"))
                # get color of selected point
                color = self._selected_idxs.get(ds).get(str(ev_idx)).get_color()
                self.axis.plot(ds.ef_grid, ef, color=color, label=label)
        self.axis.axhline(y=0, linestyle="dotted", color="grey")
        if isinstance(self.data, LegolasDataSet):
            self.axis.axvline(x=self.data.x_start, linestyle="dotted", color="grey")
        self.axis.set_title(self._get_title(ef_name))
        self.axis.legend(loc="best", fontsize=8)

    def _get_title(self, ef_name):
        """
        Creates the title of the eigenfunction plot.
        If the eigenfunction is retransformed an 'r' is prepended if appropriate,
        along with Re/Im depending on the real/imaginary part shown.
        Additionally, LaTeX formatting is used and numbers are replaced with the
        corresponding suffix: :math:`(1, 2, 3)`
        becomes :math:`(x, y, z)` or :math:`(r, \\theta, z)`.

        Parameters
        ----------
        ef_name : str
            The name of the eigenfunction as present in the datfile.

        Returns
        -------
        name : str
            The 'new' name for the eigenfunction, used as title.
        """
        part = "Re"
        if not self._use_real_part:
            part = "Im"
        suffix = ("_x", "_y", "_z")
        if self.data.geometry == "cylindrical":
            suffix = ("_r", r"_\theta", "_z")
            if ef_name in ("rho", "v1", "v3", "T", "a2") and self._retransform:
                ef_name = f"r{ef_name}"
        for i, idx in enumerate("123"):
            ef_name = ef_name.replace(idx, suffix[i])
        ef_name = ef_name.replace("rho", r"\rho")
        name = f"{part}(${ef_name}$) eigenfunction"
        return name

    def _mark_points_without_data_written(self):
        self._condition_to_make_transparent = "efs_written"
        super()._mark_points_without_data_written()
