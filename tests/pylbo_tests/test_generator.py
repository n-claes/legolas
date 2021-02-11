import pytest
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


def test_invalid_dict_type(tmp_path):
    gen = ParfileGenerator({"gridpoints": 100.5}, output_dir=tmp_path)
    with pytest.raises(TypeError):
        gen.create_namelist_from_dict()


def test_unknown_namelist_item(tmp_path):
    from pylbo.exceptions import ParfileGenerationError

    gen = ParfileGenerator({"gridpoints": 10, "unknown_item": 5}, output_dir=tmp_path)
    with pytest.raises(ParfileGenerationError):
        gen.create_namelist_from_dict()


def test_invalid_nb_runs(tmp_path):
    from pylbo.exceptions import ParfileGenerationError

    gen = ParfileGenerator(
        {"gridpoints": [50, 100], "number_of_runs": 1}, output_dir=tmp_path
    )
    with pytest.raises(ParfileGenerationError):
        gen.create_namelist_from_dict()


def test_savelist_not_present(tmp_path):
    import f90nml

    gen = ParfileGenerator(
        {"gridpoints": 100, "equilibrium_type": "adiabatic_homo"}, output_dir=tmp_path
    )
    gen.create_namelist_from_dict()
    parfile = gen.generate_parfiles()
    parfile_dict = f90nml.read(*parfile)
    assert "savelist" in parfile_dict.keys()
