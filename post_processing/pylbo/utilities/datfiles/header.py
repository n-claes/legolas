from typing import Any, BinaryIO

from pylbo._version import VersionHandler
from pylbo.visualisation.utils import ensure_attr_set


class LegolasHeader:
    """Baseclass for a Legolas header"""

    def __init__(self, istream: BinaryIO, version: VersionHandler):
        self.legolas_version = version
        self.data = {}
        self._str_len = None
        self._str_len_array = None
        self._set_str_lengths(istream)
        [ensure_attr_set(self, attr) for attr in ("_str_len", "_str_len_array")]
        self.read_header_data(istream)
        self.read_data_offsets(istream)

    def __str__(self) -> str:
        return "".join([f"{key}: {self.data.get(key)}\n" for key in self.data.keys()])

    def __getitem__(self, key: str) -> Any:
        return self.data[key]

    def _set_str_lengths(self, istream: BinaryIO) -> None:
        pass

    def get(self, key: str, default: Any = None) -> Any:
        return self.data.get(key, default)

    def read_header_data(self, istream: BinaryIO) -> None:
        raise NotImplementedError()

    def read_data_offsets(self, istream: BinaryIO) -> None:
        pass
