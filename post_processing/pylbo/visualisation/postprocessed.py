import numpy as np
from pylbo.utilities.toolbox import transform_to_numpy
from pylbo.data_containers import LegolasDataSet
from pylbo.exceptions import PostprocessedNotPresent
from pylbo.visualisation.eigenfunction_interface import EigenfunctionInterface


class PostprocessedHandler(EigenfunctionInterface):
    """
    Main handler for post-processed quantities.
    """

    def __init__(self, data, pp_ax):
        super().__init__(data, pp_ax)
        self._function_names = self.data.pp_names

    def _check_data_is_present(self):
        if not any(transform_to_numpy(self.data.pp_written)):
            raise PostprocessedNotPresent(
                "None of the given datfiles has post-processed quantities "
                "written to it."
            )

    def update_plot(self):
        self.axis.clear()
        if not self._selected_idxs:
            return
        pp_name = self._function_names[self._selected_name_idx]
        for ds, idxs_dict in self._selected_idxs.items():
            idxs = np.array([int(idx) for idx in idxs_dict.keys()])
            pp_container = ds.get_postprocessed(ev_idxs=idxs)
            for ev_idx, ppqs in zip(idxs, pp_container):
                pp = ppqs.get(pp_name)
                if self._use_real_part:
                    pp = pp.real
                else:
                    pp = pp.imag
                label = super()._get_label(ds, ev_idx, ppqs.get("eigenvalue"))
                # get color of selected point
                color = self._selected_idxs.get(ds).get(str(ev_idx)).get_color()
                self.axis.plot(ds.ef_grid, pp, color=color, label=label)
        self.axis.axhline(y=0, linestyle="dotted", color="grey")
        if isinstance(self.data, LegolasDataSet):
            self.axis.axvline(x=self.data.x_start, linestyle="dotted", color="grey")
        self.axis.set_title(self._get_title(pp_name))
        self.axis.legend(loc="best", fontsize=8)

    def _get_title(self, pp_name):
        """
        Creates the title of the post-processed quantity plot.
        LaTeX formatting is used and numbers are replaced with the
        corresponding suffix: :math:`(1, 2, 3)`
        becomes :math:`(x, y, z)` or :math:`(r, \\theta, z)`.

        Parameters
        ----------
        pp_name : str
            The name of the post-processed quantity as present in the datfile.

        Returns
        -------
        name : str
            The 'new' name for the post-processed quantity, used as title.
        """
        part = "Re"
        if not self._use_real_part:
            part = "Im"
        suffix = ("_x", "_y", "_z")
        if self.data.geometry == "cylindrical":
            suffix = ("_r", r"_\theta", "_z")
        for i, idx in enumerate("123"):
            pp_name = pp_name.replace(idx, suffix[i])
        pp_name = pp_name.replace("div", "\\nabla\\cdot")
        pp_name = pp_name.replace("curl", "\\nabla\\times")
        pp_name = pp_name.replace("para", "\\parallel")
        pp_name = pp_name.replace("perp", "\\perp")
        name = fr"{part}(${pp_name}$)"
        return name

    def _mark_points_without_data_written(self):
        self._condition_to_make_transparent = "pp_written"
        super()._mark_points_without_data_written()
