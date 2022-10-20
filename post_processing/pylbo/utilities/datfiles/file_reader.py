import struct
from functools import wraps
from typing import BinaryIO, List, Union

from pylbo._version import VersionHandler


def requires_version(version_needed, default=None):
    def check_version(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if args[0].legolas_version < version_needed:
                return default
            return func(*args, **kwargs)

        return wrapper

    return check_version


class LegolasFileReader:
    SIZE_CHAR = struct.calcsize("c")
    SIZE_INT = struct.calcsize("i")
    SIZE_BOOL = SIZE_INT  # fortran logical is 4-byte integer
    SIZE_DOUBLE = struct.calcsize("d")
    SIZE_COMPLEX = struct.calcsize(2 * "d")
    ALIGN = "="

    def __init__(self, istream: BinaryIO):
        self._istream = istream
        self.istream.seek(0)
        self.legolas_version = self._read_legolas_version()

    @property
    def istream(self) -> BinaryIO:
        """Returns the input stream of the file manager."""
        return self._istream

    def read_string_from_istream(
        self, length: int, amount: int = 1
    ) -> Union[str, List[str]]:
        """
        Reads a string from the input stream.

        Parameters
        ----------
        length : int
            The length of the string to read.
        amount : int, optional
            The amount of strings to read, by default 1

        Returns
        -------
        str, list of str
            The string or list of strings read from the input stream.
        """
        fmt = self.ALIGN + amount * length * "c"
        hdr = struct.unpack(fmt, self.istream.read(struct.calcsize(fmt)))
        if amount == 1:
            return b"".join(hdr).strip().decode()
        else:
            return [
                b"".join(hdr[i : i + length]).strip().decode()
                for i in range(0, amount * length, length)
            ]

    def _read_number_from_istream(self, kind: str, amount: int = 1) -> complex:
        """
        Reads a number from the input stream.

        Parameters
        ----------
        kind : str
            The kind of number to read.
            Can be "i" for integer, "d" for double or "c" for complex.
        amount : int, optional
            The amount of numbers to read, by default 1

        Returns
        -------
        int, float, complex
            The number or list of numbers read from the input stream.
        """
        fmt = self.ALIGN + amount * kind
        hdr = struct.unpack(fmt, self.istream.read(struct.calcsize(fmt)))
        if amount == 1:
            # for single values unpack and return single number
            (hdr,) = hdr
        return hdr

    def read_int_from_istream(self, amount: int = 1) -> int:
        """
        Reads an integer from the input stream.

        Parameters
        ----------
        amount : int, optional
            The amount of integers to read, by default 1
        """
        return self._read_number_from_istream(kind="i", amount=amount)

    def read_float_from_istream(self, amount: int = 1) -> float:
        """
        Reads a float from the input stream.

        Parameters
        ----------
        amount : int, optional
            The amount of floats to read, by default 1
        """
        return self._read_number_from_istream(kind="d", amount=amount)

    def read_complex_from_istream(self, amount: int = 1) -> complex:
        """
        Reads a complex number from the input stream.

        Parameters
        ----------
        amount : int, optional
            The amount of complex numbers to read, by default 1
        """
        return complex(*self._read_number_from_istream(kind="d", amount=2 * amount))

    def read_boolean_from_istream(self) -> bool:
        """Reads a boolean from the input stream."""
        return bool(self.read_int_from_istream())

    def _read_legolas_version(self):
        version_name = self.read_string_from_istream(length=len("legolas_version"))
        if version_name == "legolas_version":
            # formatted version is character of length 10
            version = self.read_string_from_istream(length=10)
        elif version_name == "datfile_version":
            # old numbering, single integer
            version = f"0.{str(self.read_int_from_istream())}.0"
        else:
            # very old numbering
            self.istream.seek(0)
            version = "0.0.0"
        return VersionHandler(version)
