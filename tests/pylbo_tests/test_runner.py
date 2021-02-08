import pytest
from pylbo.automation.api import run_legolas


def test_invalid_parfiles():
    with pytest.raises(FileNotFoundError):
        run_legolas(["unknown_parfile.par"])


def test_invalid_executable(tmp_path):
    # create temporary parfile
    parfile = tmp_path / "parfile.par"
    parfile.write_text("content")
    with pytest.raises(FileNotFoundError):
        run_legolas(parfile, executable="unknown_exe")
