import copy
import multiprocessing
from pathlib import Path

import pylbo
import pytest
from pylbo.automation.api import run_legolas
from pylbo.automation.runner import LegolasRunner


def test_invalid_parfiles():
    with pytest.raises(FileNotFoundError):
        run_legolas(["unknown_parfile.par"])


def test_invalid_executable(tmpdir):
    # create temporary parfile
    parfile = tmpdir / "parfile.par"
    parfile.write_text("content")
    with pytest.raises(FileNotFoundError):
        run_legolas(parfile, executable="unknown_exe")


def test_cpu_count(default_parfile):
    cpus_available = multiprocessing.cpu_count()
    run = LegolasRunner(
        parfiles=default_parfile,
        remove_parfiles=True,
        nb_cpus=cpus_available + 1,
        executable=None,
    )
    assert run.nb_cpus == cpus_available


def test_multirun(default_pf_dict):
    pf_dict = copy.deepcopy(default_pf_dict)
    pf_dict["number_of_runs"] = 4
    pf_dict["gridpoints"] = [10] * 4
    parfiles = pylbo.generate_parfiles(pf_dict)
    pylbo.run_legolas(parfiles, nb_cpus=2, remove_parfiles=True)
    for i in range(1, 5):
        filepath = (
            Path(pf_dict["output_folder"]) / f"{i:04d}{pf_dict['basename_datfile']}.dat"
        )
        assert filepath.is_file()
