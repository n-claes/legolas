import struct
from pathlib import Path
from typing import Union

import numpy as np
from pylbo.utilities.logger import pylboLogger
from pylbo.utilities.toolbox import transform_to_numpy
from pylbo.visualisation.modes.mode_data import ModeVisualisationData
from pylbo.visualisation.utils import validate_ef_name
from tqdm import tqdm


class VTKDataCube3D:
    """
    Class to prepare eigenmode solution data for export to VTK files.

    Parameters
    ----------
    data : ModeVisualisationData
        The data for the visualisation.
    u2 : np.ndarray
        The :math:`y` or :math:`\\theta` coordinates of the eigenmode solution.
    u3 : np.ndarray
        The :math:`z` coordinates of the eigenmode solution.

    Attributes
    ----------
    data : ModeVisualisationData
        The data for the visualisation.
    dims : tuple
        The dimensions of the eigenmode solution.
    u1_data : np.ndarray
        The :math:`x` coordinates of the eigenmode solution.
    u2_data : np.ndarray
        The :math:`y` or :math:`\\theta` coordinates of the eigenmode solution.
    u3_data : np.ndarray
        The :math:`z` coordinates of the eigenmode solution.
    _vtk_dtype : str
        The VTK data type, defaults to "float".
    _vtk_byte_order : str
        The VTK byte order, defaults to ">" (big endian).
    _vtk_fmt : str
        The VTK data format, defaults to ">f".
    _datacube : np.ndarray
        The eigenmode solution data in a 3D datacube ordered as :math:`[x, y, z]` in
        Cartesian coordinates or :math:`[r, \\theta, z]` in cylindrical coordinates.
    """

    def __init__(
        self,
        data: ModeVisualisationData,
        u2: np.ndarray,
        u3: np.ndarray,
    ) -> None:
        self.data = data
        self._u1 = data.ds.ef_grid
        self._u2 = u2
        self._u3 = u3
        self._vtk_dtype = None
        self._vtk_byte_order = ">"  # big endian
        self._vtk_fmt = None
        self._datacube = None
        self._pbar = None

        self.dims = tuple([len(self._u1), len(self._u2), len(self._u3)])
        geometry = self.data.ds.geometry.lower()
        pylboLogger.info(f"creating {geometry} datacube with dimensions {self.dims}")
        if geometry == "cartesian":
            self.u1_data, self.u2_data, self.u3_data = np.meshgrid(
                self._u1, self._u2, self._u3, indexing="ij"
            )
        else:
            r_3d, theta_3d, self.u3_data = np.meshgrid(
                self._u1, self._u2, self._u3, indexing="ij"
            )
            self.u1_data = r_3d * np.cos(theta_3d)
            self.u2_data = r_3d * np.sin(theta_3d)

    def _validate_and_set_dtype(self, dtype: str) -> None:
        """
        Validates and sets the VTK data type.

        Parameters
        ----------
        dtype : str
            The VTK data type. Valid values are "float32" and "float64".

        Raises
        ------
        ValueError
            If the VTK data type is not valid.
        """
        if dtype == "float32":
            self._vtk_dtype = "float"
            self._vtk_fmt = f"{self._vtk_byte_order}f"
        elif dtype == "float64":
            self._vtk_dtype = "double"
            self._vtk_fmt = f"{self._vtk_byte_order}d"
        else:
            raise ValueError(f"dtype {dtype} not supported.")

    def _log_info(self, msg: str) -> None:
        if self._pbar is not None:
            return
        pylboLogger.info(msg)

    def _update_pbar(self) -> None:
        if self._pbar is not None:
            self._pbar.update()

    def get_solution(self, ef_name: str, time: float) -> np.ndarray:
        """
        Parameters
        ----------
        ef_name : str
            The name of the eigenfunction to use.
        time : float
            The time at which to return the eigenmode solution.

        Returns
        -------
        np.ndarray
            The eigenmode solution at the given time.
        """
        name = validate_ef_name(self.data.ds, ef_name)
        ef_3d = np.broadcast_to(
            self.data._efs.get(name),
            shape=(len(self._u3), len(self._u2), len(self._u1)),
        ).transpose()
        return self.data.get_mode_solution(ef_3d, self.u2_data, self.u3_data, time)

    def export_to_vtk(
        self,
        filename: str,
        time: Union[float, np.ndarray],
        names: Union[str, list[str]],
        dtype: str = "float32",
    ) -> None:
        """
        Exports the current datacube to a VTK file at the given time(s).

        Parameters
        ----------
        filename : str
            The name of the VTK file to export to, no extension needed.
        time : float
            The time at which to export the eigenmode solution.
        names : str or list[str]
            The name(s) of the eigenfunction(s) to export. If a list is given,
            the eigenmode solutions will be exported as multiple scalar fields in
            the VTK file.
        dtype : str, optional
            The VTK data type, defaults to "float32" (32 bit floating point).
            Can be set to "float64" (64 bit floating point) but uses more memory.
        """
        time = transform_to_numpy(time)
        names = transform_to_numpy(names)
        self._validate_and_set_dtype(dtype)
        filename = Path(filename).with_suffix("")  # remove extension
        if len(time) > 1:
            self._pbar = tqdm(total=len(time), desc="writing VTK files", unit="file")
        for it, t in enumerate(time):
            vtkfile = Path(f"{filename}_t{it:04d}.vtk")
            self._log_info(f"writing {vtkfile.name}...")
            self._write_vtk_file(vtkfile, t, names)
            self._log_info(f"export to VTK completed, stored at {vtkfile.resolve()}")
            self._update_pbar()

    def _write_vtk_file(self, vtkfile: str, time: float, names: np.ndarray) -> None:
        """
        Writes the VTK file.

        Parameters
        ----------
        vtkfile : str
            The name of the VTK file to write to.
        time : float
            The time at which to export the eigenmode solution.
        names : np.ndarray
            The name(s) of the eigenfunction(s) to write.
        """
        self._log_info("writing VTK file header...")
        self._write_vtk_header(vtkfile)
        self._log_info("writing VTK structured grid data...")
        self._write_vtk_grid_data(vtkfile)
        self._log_info("writing VTK point data...")
        self._write_vtk_point_data(vtkfile, time, names)

    def _write_vtk_header(self, vtkfile):
        """
        Writes the VTK file header. This includes the VTK file version, the
        data type, the grid dimensions and the number of points.

        Parameters
        ----------
        vtkfile : str
            The name of the VTK file to write to.
        """
        with open(vtkfile, "w") as ostream:
            vtktitle = vtkfile.stem if len(vtkfile.stem) < 256 else "vtk output"
            ostream.write("# vtk DataFile Version 3.0 \n")
            ostream.write(f"{vtktitle} \n")
            ostream.write("BINARY \n")
            ostream.write("DATASET STRUCTURED_GRID \n")
            ostream.write(f"DIMENSIONS {self.dims[0]} {self.dims[1]} {self.dims[2]} \n")
            ostream.write(f"POINTS {np.prod(self.dims)} {self._vtk_dtype} \n")

    def _write_vtk_grid_data(self, vtkfile):
        """
        Writes the VTK grid coordinates.

        Parameters
        ----------
        vtkfile : str
            The name of the VTK file to write to.
        """
        with open(vtkfile, "ab") as ostream:
            for k in range(self.dims[2]):
                for j in range(self.dims[1]):
                    for i in range(self.dims[0]):
                        ostream.write(struct.pack(self._vtk_fmt, self.u1_data[i, j, k]))
                        ostream.write(struct.pack(self._vtk_fmt, self.u2_data[i, j, k]))
                        ostream.write(struct.pack(self._vtk_fmt, self.u3_data[i, j, k]))

    def _write_vtk_point_data(self, vtkfile, time, names):
        """
        Writes the VTK point data as a series of scalar fields, one for each
        eigenfunction name given.

        Parameters
        ----------
        vtkfile : str
            The name of the VTK file to write to.
        time : float
            The time at which to export the eigenmode solution.
        names : np.ndarray
            The name(s) of the eigenfunction(s) to write.
        """
        with open(vtkfile, "a") as ostream:
            ostream.write(f"\nPOINT_DATA {np.prod(self.dims)} \n")
        for name in names:
            with open(vtkfile, "a") as ostream:
                ostream.write(f"SCALARS {name} {self._vtk_dtype} \n")
                ostream.write("LOOKUP_TABLE default \n")
            solution = self.get_solution(name, time)
            with open(vtkfile, "ab") as ostream:
                for k in range(self.dims[2]):
                    for j in range(self.dims[1]):
                        for i in range(self.dims[0]):
                            ostream.write(struct.pack(self._vtk_fmt, solution[i, j, k]))
