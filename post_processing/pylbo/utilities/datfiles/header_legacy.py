from __future__ import annotations

from typing import Any

import numpy as np
from pylbo.utilities.datfiles.file_manager import LegolasFileManager, requires_version


class LegolasLegacyHeader:
    def __init__(self, filemanager: LegolasFileManager) -> None:
        self.legolas_version = filemanager.legolas_version
        self.manager = filemanager
        self.data = {}

        self._str_len = filemanager.read_int_from_istream()
        self._str_len_array = filemanager.read_int_from_istream()
        self.read_header_data()
        self.read_data_offsets()

    def __str__(self) -> str:
        return "".join([f"{key}: {self.data.get(key)}\n" for key in self.data.keys()])

    def __getitem__(self, key: str) -> Any:
        return self.data[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.data.get(key, default)

    def read_header_data(self) -> None:
        data = {}
        manager = self.manager

        data["geometry"] = manager.read_string_from_istream(length=self._str_len)
        data["x_start"], data["x_end"] = manager.read_float_from_istream(amount=2)

        for key in ("", "gauss_", "matrix_", "ef_"):
            data[f"{key}gridpoints"] = manager.read_int_from_istream()

        data["gamma"] = manager.read_float_from_istream()
        data["eq_type"] = manager.read_string_from_istream(length=self._str_len)

        data["has_efs"] = manager.read_boolean_from_istream()
        data["has_derived_efs"] = self._read_has_derived_efs()
        data["has_matrices"] = manager.read_boolean_from_istream()
        data["has_eigenvectors"] = self._read_has_eigenvectors()
        data["has_residuals"] = self._read_has_residuals()
        (
            data["ef_subset_used"],
            data["ef_subset_center"],
            data["ef_subset_radius"],
        ) = self._read_ef_subset_properties()
        data["parameters"] = self._read_parameters()
        data["equilibrium_names"] = self._read_equilibrium_names()

        data["units"] = self._read_units()
        data["nb_eigenvalues"] = (
            manager.read_int_from_istream()
            if self.legolas_version >= "1.0.2"
            else data["matrix_gridpoints"]
        )
        data["offsets"] = {}
        self.data.update(data)

    def read_data_offsets(self) -> None:
        istream = self.manager.istream

        offsets = {}
        # eigenvalue offset
        offsets["eigenvalues"] = istream.tell()
        bytesize = self.data["nb_eigenvalues"] * self.manager.SIZE_COMPLEX
        istream.seek(istream.tell() + bytesize)
        # grid offset
        offsets["grid"] = istream.tell()
        bytesize = self.data["gridpoints"] * self.manager.SIZE_DOUBLE
        istream.seek(istream.tell() + bytesize)
        # grid gauss offset
        offsets["grid_gauss"] = istream.tell()
        bytesize = self.data["gauss_gridpoints"] * self.manager.SIZE_DOUBLE
        istream.seek(istream.tell() + bytesize)
        # equilibrium arrays offset
        offsets["equilibrium_arrays"] = istream.tell()
        bytesize = (
            self.data["gauss_gridpoints"]
            * len(self.data["equilibrium_names"])
            * self.manager.SIZE_DOUBLE
        )
        istream.seek(istream.tell() + bytesize)

        offsets.update(self._get_eigenfunction_offsets())
        offsets.update(self._get_derived_eigenfunction_offsets())
        offsets.update(self._get_eigenvector_offsets())
        offsets.update(self._get_residuals_offsets())
        offsets.update(self._get_matrices_offsets())

        self.data["offsets"].update(offsets)

    @requires_version("1.1.3", default=False)
    def _read_has_derived_efs(self) -> bool:
        return self.manager.read_boolean_from_istream()

    @requires_version("1.3.0", default=False)
    def _read_has_eigenvectors(self) -> bool:
        return self.manager.read_boolean_from_istream()

    @requires_version("1.3.0", default=False)
    def _read_has_residuals(self) -> bool:
        return self.manager.read_boolean_from_istream()

    @requires_version("1.1.4", default=(False, None, None))
    def _read_ef_subset_properties(self) -> tuple(bool, complex, float):
        used = self.manager.read_boolean_from_istream()
        center = self.manager.read_complex_from_istream()
        radius = self.manager.read_float_from_istream()
        return (used, center, radius)

    def _read_parameters(self) -> dict:
        nb_params = self.manager.read_int_from_istream()
        len_param_name = (
            self.manager.read_int_from_istream()
            if self.legolas_version >= "1.0.2"
            else self._str_len_array
        )
        parameter_names = self.manager.read_string_from_istream(
            length=len_param_name, amount=nb_params
        )
        parameter_values = self.manager.read_float_from_istream(amount=nb_params)
        return {
            name: value
            for name, value in zip(parameter_names, parameter_values)
            if not np.isnan(value)
        }

    def _read_equilibrium_names(self) -> list[str]:
        nb_names = self.manager.read_int_from_istream()
        len_name = (
            self.manager.read_int_from_istream()
            if self.legolas_version >= "1.0.2"
            else self._str_len_array
        )
        return self.manager.read_string_from_istream(length=len_name, amount=nb_names)

    def _read_units(self) -> dict:
        units = {"cgs": self.manager.read_boolean_from_istream()}
        if self.legolas_version >= "1.0.2":
            nb_units, len_unit_name = self.manager.read_int_from_istream(amount=2)
            unit_names = self.manager.read_string_from_istream(
                length=len_unit_name, amount=nb_units
            )
        else:
            unit_names = [
                "unit_length",
                "unit_time",
                "unit_density",
                "unit_velocity",
                "unit_temperature",
                "unit_pressure",
                "unit_magneticfield",
                "unit_numberdensity",
                "unit_lambdaT",
                "unit_conduction",
                "unit_resistivity",
            ]
            nb_units = len(unit_names)
        unit_values = self.manager.read_float_from_istream(amount=nb_units)
        for name, value in zip(unit_names, unit_values):
            units[name] = value
        # mean molecular weight is added in 1.1.2, before this it defaults to 1
        units.setdefault("mean_molecular_weight", 1.0)
        return units

    def _get_eigenfunction_offsets(self) -> dict:
        if not self.data["has_efs"]:
            return {}
        # eigenfunction names
        nb_efs = self.manager.read_int_from_istream()
        self.data["ef_names"] = self.manager.read_string_from_istream(
            length=self._str_len_array, amount=nb_efs
        )
        # eigenfunction grid offset
        offsets = {"ef_grid": self.manager.istream.tell()}
        bytesize = self.data["ef_gridpoints"] * self.manager.SIZE_DOUBLE
        self.manager.istream.seek(self.manager.istream.tell() + bytesize)
        # ef written flags
        self._set_ef_written_flags()
        # eigenfunction offsets
        offsets["ef_arrays"] = self.manager.istream.tell()
        bytesize = (
            self.data["ef_gridpoints"]
            * len(self.data["ef_written_idxs"])
            * nb_efs
            * self.manager.SIZE_COMPLEX
        )
        self.manager.istream.seek(self.manager.istream.tell() + bytesize)
        return offsets

    def _set_ef_written_flags(self) -> None:
        if self.legolas_version < "1.1.4":
            self.data["ef_written_flags"] = np.asarray(
                [True] * self.data["nb_eigenvalues"], dtype=bool
            )
            self.data["ef_written_idxs"] = np.arange(0, self.data["nb_eigenvalues"])
            return

        ef_flags_size = self.manager.read_int_from_istream()
        self.data["ef_written_flags"] = np.asarray(
            self.manager.read_int_from_istream(amount=ef_flags_size), dtype=bool
        )
        ef_idxs_size = self.manager.read_int_from_istream()
        self.data["ef_written_idxs"] = (
            np.asarray(
                self.manager.read_int_from_istream(amount=ef_idxs_size), dtype=int
            )
            - 1
        )  # -1 to correct for Fortran 1-based indexing
        # sanity check
        assert all(
            self.data["ef_written_idxs"] == np.where(self.data["ef_written_flags"])[0]
        )

    def _get_derived_eigenfunction_offsets(self) -> dict:
        if not self.data["has_derived_efs"]:
            return {}
        nb_defs = self.manager.read_int_from_istream()
        self.data["derived_ef_names"] = self.manager.read_string_from_istream(
            length=self._str_len_array, amount=nb_defs
        )
        offsets = {"derived_ef_arrays": self.manager.istream.tell()}
        bytesize = (
            self.data["ef_gridpoints"]
            * len(self.data["ef_written_idxs"])
            * nb_defs
            * self.manager.SIZE_COMPLEX
        )
        self.manager.istream.seek(self.manager.istream.tell() + bytesize)
        return offsets

    def _get_eigenvector_offsets(self) -> dict:
        if not self.data["has_eigenvectors"]:
            return {}
        len_eigvecs, nb_eigvecs = self.manager.read_int_from_istream(amount=2)
        offsets = {
            "eigenvectors": self.manager.istream.tell(),
            "eigenvector_length": len_eigvecs,
            "nb_eigenvectors": nb_eigvecs,
        }
        bytesize = len_eigvecs * nb_eigvecs * self.manager.SIZE_COMPLEX
        self.manager.istream.seek(self.manager.istream.tell() + bytesize)
        return offsets

    def _get_residuals_offsets(self) -> dict:
        if not self.data["has_residuals"]:
            return {}
        nb_residuals = self.manager.read_int_from_istream()
        offsets = {
            "residuals": self.manager.istream.tell(),
            "nb_residuals": nb_residuals,
        }
        bytesize = nb_residuals * self.manager.SIZE_DOUBLE
        self.manager.istream.seek(self.manager.istream.tell() + bytesize)
        return offsets

    def _get_matrices_offsets(self) -> dict:
        if not self.data["has_matrices"]:
            return {}
        nonzero_B_elements = self.manager.read_int_from_istream()
        nonzero_A_elements = self.manager.read_int_from_istream()
        # B offsets, written as (row, column, value)
        byte_size = (
            2 * self.manager.SIZE_INT + self.manager.SIZE_DOUBLE
        ) * nonzero_B_elements
        offsets = {
            "matrix_B": self.manager.istream.tell(),
            "nonzero_B_elements": nonzero_B_elements,
        }
        self.manager.istream.seek(self.manager.istream.tell() + byte_size)
        # A offsets
        offsets["matrix_A"] = self.manager.istream.tell()
        offsets["nonzero_A_elements"] = nonzero_A_elements
        return offsets
