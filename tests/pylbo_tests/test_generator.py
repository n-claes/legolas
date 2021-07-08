import f90nml
import pytest
import numpy as np
from pathlib import Path
from pylbo.automation.generator import ParfileGenerator


def test_basename_none(tmp_path):
    gen = ParfileGenerator({}, output_dir=tmp_path)
    assert gen.basename == "parfile"


def test_outputdir_none():
    parfiledir = Path.cwd() / "parfiles"
    gen = ParfileGenerator({})
    assert gen.output_dir == parfiledir
    parfiledir.rmdir()


def test_outputdir_invalid():
    with pytest.raises(NotADirectoryError):
        ParfileGenerator({}, output_dir="invalid_dir")


def test_invalid_dict_type():
    gen = ParfileGenerator({"gridpoints": 100.5})
    with pytest.raises(TypeError):
        gen.create_namelist_from_dict()


def test_unknown_namelist_item():
    from pylbo.exceptions import ParfileGenerationError

    gen = ParfileGenerator({"gridpoints": 10, "unknown_item": 5})
    with pytest.raises(ParfileGenerationError):
        gen.create_namelist_from_dict()


def test_invalid_nb_runs():
    from pylbo.exceptions import ParfileGenerationError

    gen = ParfileGenerator({"gridpoints": [50, 100], "number_of_runs": 1})
    with pytest.raises(ParfileGenerationError):
        gen.create_namelist_from_dict()


def test_savelist_not_present(tmp_path):
    gen = ParfileGenerator(
        {"gridpoints": 100, "equilibrium_type": "adiabatic_homo"}, output_dir=tmp_path
    )
    gen.create_namelist_from_dict()
    parfile = gen.generate_parfiles()
    parfile_dict = f90nml.read(*parfile)
    assert "savelist" in parfile_dict.keys()


def test_paramlist_present():
    gen = ParfileGenerator(
        {
            "gridpoints": 50,
            "equilibrium_type": "adiabatic_homo",
            "parameters": {"k2": 2, "k3": 0, "cte_T0": 5},
        },
    )
    gen.create_namelist_from_dict()
    assert "paramlist" in gen.container.keys()
    for key, value in zip(("k2", "k3", "cte_T0"), ([2], [0], [5])):
        assert gen.container["paramlist"].get(key) == value


def test_nb_runs():
    gen = ParfileGenerator(
        {
            "gridpoints": [50, 80, 150],
            "number_of_runs": 3,
            "equilibrium_type": "adiabatic_homo",
        },
    )
    gen.create_namelist_from_dict()
    assert gen.container["gridlist"].get("gridpoints") == [50, 80, 150]
    assert (
        gen.container["equilibriumlist"].get("equilibrium_type")
        == ["adiabatic_homo"] * 3
    )


def test_sigma_complex():
    gen = ParfileGenerator({"solver": "arnoldi", "sigma": 5.0})
    gen.create_namelist_from_dict()
    sigma = gen.container["solvelist"].get("sigma")[0]
    assert isinstance(sigma, complex)
    assert np.isclose(sigma, 5.0 + 0j)


def test_logfile_name(tmp_path):
    gen = ParfileGenerator(
        {
            "gridpoints": 50,
            "equilibrium_type": "user_defined",
            "basename_logfile": "mylog",
        },
        output_dir=tmp_path,
    )
    gen.create_namelist_from_dict()
    parfile = gen.generate_parfiles()
    parfile_dict = f90nml.read(*parfile)
    assert parfile_dict["savelist"].get("basename_logfile") == "mylog"
