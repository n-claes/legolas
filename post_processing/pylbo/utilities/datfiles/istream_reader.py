from __future__ import annotations

import struct
from functools import wraps
from typing import BinaryIO, Union

DTYPES = {
    "str": "c",
    "int": "i",
    "bool": "i",
    "float": "d",
    "complex": 2 * "d",
}
NP_DTYPES = {
    "str": str,
    "int": int,
    "bool": bool,
    "float": float,
    "complex": complex,
}
BYTE_ORDERS = {
    "little": "<",
    "big": ">",
    "native": "=",
}

SIZE_CHAR = struct.calcsize(DTYPES["str"])
SIZE_INT = struct.calcsize(DTYPES["int"])
SIZE_BOOL = struct.calcsize(DTYPES["bool"])
SIZE_DOUBLE = struct.calcsize(DTYPES["float"])
SIZE_COMPLEX = struct.calcsize(DTYPES["complex"])


def requires_version(version_needed, default=None):
    def check_version(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if args[0].legolas_version < version_needed:
                return default
            return func(*args, **kwargs)

        return wrapper

    return check_version


def read_string_from_istream(
    istream: BinaryIO,
    length: int,
    amount: int = 1,
    offset: int = None,
    byte_order: str = "native",
) -> Union[str, list[str]]:
    """
    Reads a string from the input stream.

    Parameters
    ----------
    istream : BinaryIO
        The input stream to read from.
    length : int
        The length of the string to read.
    amount : int, optional
        The amount of strings to read, by default 1.
    offset : int, optional
        The offset to seek to before reading, by default `None`.
    byte_order : str, optional
        The byte order to use, by default "native".

    Returns
    -------
    str, list of str
        The string(s) read from the input stream.
    """
    if offset is not None:
        istream.seek(offset)
    fmt = BYTE_ORDERS[byte_order] + amount * length * DTYPES["str"]
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    if amount == 1:
        return b"".join(hdr).strip().decode()
    return [
        b"".join(hdr[i : i + length]).strip().decode()
        for i in range(0, amount * length, length)
    ]


def read_int_from_istream(
    istream: BinaryIO,
    amount: int = 1,
    offset: int = None,
    byte_order: str = "native",
) -> Union[int, tuple[int, ...]]:
    """
    Reads an integer from the input stream.

    Parameters
    ----------
    istream : BinaryIO
        The input stream to read from.
    amount : int, optional
        The amount of integers to read, by default 1.
    offset : int, optional
        The offset to seek to before reading, by default `None`.
    byte_order : str, optional
        The byte order to use, by default "native".

    Returns
    -------
    int, tuple of int
        The integer(s) read from the input stream.
    """
    if offset is not None:
        istream.seek(offset)
    fmt = BYTE_ORDERS[byte_order] + amount * DTYPES["int"]
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    if amount == 1:
        (hdr,) = hdr  # unpack for single values
    return hdr


def read_boolean_from_istream(
    istream: BinaryIO,
    offset: int = None,
    byte_order: str = "native",
) -> bool:
    """
    Reads a boolean from the input stream.

    Parameters
    ----------
    istream : BinaryIO
        The input stream to read from.
    offset : int, optional
        The offset to seek to before reading, by default `None`.
    byte_order : str, optional
        The byte order to use, by default "native".

    Returns
    -------
    bool
        The boolean read from the input stream.
    """
    return bool(read_int_from_istream(istream, offset=offset, byte_order=byte_order))


def read_float_from_istream(
    istream: BinaryIO,
    amount: int = 1,
    offset: int = None,
    byte_order: str = "native",
) -> Union[float, tuple[float, ...]]:
    """
    Reads a float from the input stream.

    Parameters
    ----------
    istream : BinaryIO
        The input stream to read from.
    amount : int, optional
        The amount of floats to read, by default 1.
    offset : int, optional
        The offset to seek to before reading, by default `None`.
    byte_order : str, optional
        The byte order to use, by default "native".

    Returns
    -------
    float, tuple of float
        The float(s) read from the input stream.
    """
    if offset is not None:
        istream.seek(offset)
    fmt = BYTE_ORDERS[byte_order] + amount * DTYPES["float"]
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    if amount == 1:
        (hdr,) = hdr  # unpack for single values
    return hdr


def read_complex_from_istream(
    istream: BinaryIO,
    amount: int = 1,
    offset: int = None,
    byte_order: str = "native",
) -> Union[complex, tuple[complex, ...]]:
    """
    Reads a complex from the input stream.

    Parameters
    ----------
    istream : BinaryIO
        The input stream to read from.
    amount : int, optional
        The amount of complex numbers to read, by default 1.
    offset : int, optional
        The offset to seek to before reading, by default `None`.
    byte_order : str, optional
        The byte order to use, by default "native".

    Returns
    -------
    complex, tuple of complex
        The complex number(s) read from the input stream.
    """
    if offset is not None:
        istream.seek(offset)
    fmt = BYTE_ORDERS[byte_order] + amount * DTYPES["complex"]
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    if amount == 1:
        return complex(*hdr)  # unpack for single values
    reals, imags = hdr[::2], hdr[1::2]
    return tuple([complex(x, y) for x, y in zip(reals, imags)])


def read_mixed_from_istream(
    istream: BinaryIO,
    fmt: str,
    amount: int = 1,
    offset: int = None,
    byte_order: str = "native",
) -> tuple(complex, ...):
    """
    Reads a number of mixed types from the input stream.

    Parameters
    ----------
    istream : BinaryIO
        The input stream to read from.
    fmt : str
        The format string to use.
    amount : int, optional
        The amount of mixed types to read, by default 1.
    offset : int, optional
        The offset to seek to before reading, by default `None`.
    byte_order : str, optional
        The byte order to use, by default "native".

    Returns
    -------
    tuple of mixed
        The mixed types read from the input stream.
    """
    for char in fmt:
        if char not in DTYPES.values():
            raise ValueError(
                f"Invalid format character {char}, expected one of {DTYPES.values()}"
            )
    if offset is not None:
        istream.seek(offset)
    fmt = BYTE_ORDERS[byte_order] + amount * fmt
    return struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
