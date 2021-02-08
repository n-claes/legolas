import pytest
from pylbo._version import VersionHandler


def test_invalid_version():
    with pytest.raises(ValueError):
        VersionHandler(1.5)


def test_version_smaller():
    assert VersionHandler("1.2") < "1.3"
    assert VersionHandler("1.1.3") <= "1.1.3" <= "1.2.4"


def test_version_larger():
    assert VersionHandler("1.8") > "1.4"
    assert VersionHandler("2.1") >= "2.1" >= "2.0.5"


def test_version_equal():
    assert VersionHandler("1.2") == "1.2"


def test_version_notequal():
    assert VersionHandler("1.3") != "1.2"
