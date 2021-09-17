import shutil
import os
import subprocess
import yaml
from pathlib import Path
from argparse import ArgumentParser
from tempfile import mkstemp

DOCS = Path(__file__).resolve().parent
SRC = (Path(__file__).resolve().parents[1] / "src").resolve()
PYLBO = (Path(__file__).resolve().parents[1] / "post_processing/pylbo").resolve()
assert DOCS.is_dir()
assert SRC.is_dir()
assert PYLBO.is_dir()
FORD_FILENAME = "build_ford.md"


def _read_legolas_version():
    version_filepath = SRC / "mod_version.f08"
    with open(version_filepath) as f:
        for line in f.readlines():
            if "LEGOLAS_VERSION" in line.strip():
                return line.split("=")[-1].strip().strip('"')
    return None


def _read_pylbo_version():
    version_filepath = PYLBO / "_version.py"
    with open(version_filepath) as f:
        for line in f.readlines():
            if line.strip().startswith("__version__"):
                return line.split("=")[-1].strip().strip('"')
    return None


def parse_command_arguments():
    parser = ArgumentParser()
    parser.add_argument("branch", type=str, nargs="?", default=None)
    parser.add_argument("--clean", action="store_true")
    args = parser.parse_args()
    branch = args.branch
    clean_files = args.clean
    if branch is not None and branch not in ("stable", "develop"):
        raise ValueError(f"expected 'stable' or 'develop' as argument, got '{branch}'")
    return branch, clean_files


def generate_ford_docs(branch):
    ford_config = (DOCS / FORD_FILENAME).resolve()
    filehandle, abspath = mkstemp()
    with os.fdopen(filehandle, "w") as new_file:
        with open(ford_config) as original_file:
            for line in original_file:
                if "project:" in line:
                    new_file.write(
                        f"project: Legolas v. {_read_legolas_version()} - {branch}\n"
                    )
                else:
                    new_file.write(line)
    # copy file permissions
    shutil.copymode(ford_config, abspath)
    # remove old file
    os.remove(ford_config)
    # save new file
    shutil.move(abspath, ford_config)
    # generate documentation
    subprocess.check_call(["ford", FORD_FILENAME])
    # restore original file
    subprocess.check_call(["git", "checkout", FORD_FILENAME])


def generate_sphinx_docs(branch):
    # create temporary pylbo version file to read in with conf.py
    pylbo_tmp = (DOCS / "pylbo_version.txt").resolve()
    if pylbo_tmp.is_file():
        os.remove(pylbo_tmp)
    with open(pylbo_tmp, "w") as ostream:
        ostream.write(f"version: {_read_pylbo_version()}\n")
        ostream.write(f"release: {branch}")
    # generate documentation
    subprocess.check_call(["sphinx-build", "-b", "html", "sphinx_source", "sphinx"])
    # remove temporary file
    os.remove(pylbo_tmp)


def modify_config_yml(branch):
    if branch is None:
        return
    config_yml = (DOCS / "_config.yml").resolve()
    assert config_yml.is_file()
    with open(config_yml, "r") as file:
        config = yaml.safe_load(file)
    config["subtitle"] = f"v. {_read_legolas_version()} - {branch}"
    with open(config_yml, "w") as file:
        yaml.safe_dump(config, file, default_flow_style=False, sort_keys=False)


def modify_navigation_yml(branch):
    if branch is None:
        return
    nav_yml = (DOCS / "_data/navigation.yml").resolve()
    assert nav_yml.is_file()
    with open(nav_yml, "r") as file:
        nav = yaml.safe_load(file)
    if branch == "stable":
        button = {"title": "goto: dev ▶", "url": "https://dev.legolas.science"}
    else:
        button = {"title": "◀ goto: stable", "url": "https://legolas.science"}
    nav["main"].append(button)
    with open(nav_yml, "w") as file:
        yaml.safe_dump(nav, file, default_flow_style=False, sort_keys=False)


def main():
    branch, clean_files = parse_command_arguments()
    if branch is not None:
        branchdir = (DOCS / branch).resolve()
        if branchdir.is_dir():
            shutil.rmtree(branchdir)
        branchdir.mkdir()

    generate_ford_docs(branch)
    generate_sphinx_docs(branch)
    modify_navigation_yml(branch)
    modify_config_yml(branch)


if __name__ == "__main__":
    main()
