import numpy as np
from pylbo.utilities.toolbox import transform_to_numpy
from pylbo.data_containers import LegolasDataSet
from pylbo.exceptions import EigenfunctionsNotPresent
from pylbo.visualisation.eigenfunction_interface import EigenfunctionInterface


class DerivedEigenfunctionHandler(EigenfunctionInterface):
    """
    Main handler for derived eigenfunctions.
    """

    def __init__(self, data, def_ax, spec_ax):
        super().__init__(data, def_ax, spec_ax)
        self._function_names = self.data.derived_ef_names

    def _check_data_is_present(self):
        if not any(transform_to_numpy(self.data.derived_efs_written)):
            raise EigenfunctionsNotPresent(
                "None of the given datfiles has derived eigenfunctions "
                "written to it."
            )

    def update_plot(self):
        self.axis.clear()
        if not self._selected_idxs:
            return
        ef_name = self._function_names[self._selected_name_idx]
        for ds, idxs_dict in self._selected_idxs.items():
            idxs = np.array([int(idx) for idx in idxs_dict.keys()])
            ef_container = ds.get_derived_eigenfunctions(ev_idxs=idxs)
            for ev_idx, efs in zip(idxs, ef_container):
                ef = efs.get(ef_name)
                if self._use_real_part:
                    ef = ef.real
                else:
                    ef = ef.imag
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
        Creates the title of the derived eigenfunction plot.
        LaTeX formatting is used and numbers are replaced with the
        corresponding suffix: :math:`(1, 2, 3)`
        becomes :math:`(x, y, z)` or :math:`(r, \\theta, z)`.

        Parameters
        ----------
        ef_name : str
            The name of the derived eigenfunction as present in the datfile.

        Returns
        -------
        name : str
            The 'new' name for the derived eigenfunction, used as title.
        """
        part = "Re"
        if not self._use_real_part:
            part = "Im"
        suffix = ("_x", "_y", "_z")
        if self.data.geometry == "cylindrical":
            suffix = ("_r", r"_\theta", "_z")
        for i, idx in enumerate("123"):
            ef_name = ef_name.replace(idx, suffix[i])
        ef_name = ef_name.replace("div", "\\nabla\\cdot")
        ef_name = ef_name.replace("curl", "\\nabla\\times")
        ef_name = ef_name.replace("para", "\\parallel")
        ef_name = ef_name.replace("perp", "\\perp")
        name = fr"{part}(${ef_name}$)"
        return name

    def _mark_points_without_data_written(self):
        self._condition_to_make_transparent = "derived_efs_written"
        super()._mark_points_without_data_written()
