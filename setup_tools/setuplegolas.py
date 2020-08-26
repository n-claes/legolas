#!/usr/bin/python3
import os
import sys
import shutil
from tempfile import mkstemp
from pathlib import Path

try:
    LEGOLASDIR = os.environ['LEGOLASDIR']
except KeyError as e:
    print("The environment variable $LEGOLASDIR has not been set! Set it with: \n"
          "export $LEGOLASDIR=path_to_legolas")
    sys.exit(1)

LEGOLASDIR = Path(LEGOLASDIR).resolve()


def copy_file(filename, msg=None, location=None):
    if msg is None:
        msg = filename
    if Path(filename).resolve().is_file():
        return
    print(f">> Copying {msg} to present directory.")
    if location is None:
        srcloc = LEGOLASDIR
    else:
        srcloc = (LEGOLASDIR / location).resolve()
    path_to_file = (srcloc / filename).resolve()
    path_to_dest = (Path(os.getcwd()) / filename).resolve()
    shutil.copyfile(path_to_file, path_to_dest)


def ask_to_copy_file(filename, msg, location=None):
    if Path(filename).resolve().exists():
        return
    answer = input(f"{msg} [y/n] ").lower()
    if answer in ('yes', 'y'):
        copy_file(filename, location)


def main():
    # copy over CMake file
    copy_file("CMakeLists.txt")
    # ask to copy quick-build file
    ask_to_copy_file(filename="buildlegolas.sh",
                     msg=">> Copy over quick-build shell script?")
    # ask to copy wrapper
    ask_to_copy_file(filename="pylbo_wrapper.py",
                     msg=">> You'll need the Python wrapper to plot results after running. Copy over?")

    # check parfile
    if not list(Path(os.getcwd()).glob("*.par")):
        ask_to_copy_file("legolas_config.par",
                         msg=">> No parfile found. Copy default template?")

    # check user-defined submodule
    usr_mod_present = True
    usr_mod_filepath = (Path(os.getcwd()) / "smod_user_defined.f08").resolve()
    if not usr_mod_filepath.exists():
        usr_mod_present = False
        answer = input(">> No user-defined submodule found, copy default template? [y/n] ").lower()
        if answer in ('yes', 'y'):
            copy_file("smod_user_defined.f08", location='src/equilibria')
            usr_mod_present = True
        else:
            print(">> Submodule not copied, using default in legolas src directory.")
            # at this point we're done, and won't use the submodule
            exit()

    # if user-defined submodule is present, we have to tell CMake to use THIS one
    cmakefile = "CMakeLists.txt"
    file_handler, abspath = mkstemp()
    with os.fdopen(file_handler, "w") as new_file:
        with open(cmakefile) as original_file:
            for line in original_file:
                expr = "set(USR_SMOD_LOC equilibria)"
                if expr in line:
                    new_line = "".join([expr.replace("equilibria", f"{os.getcwd()}"), "\n"])
                else:
                    new_line = line
                new_file.write(new_line)
    # copy permissions
    shutil.copymode(cmakefile, abspath)
    # remove old file
    os.remove(cmakefile)
    # save new file
    shutil.move(abspath, cmakefile)


if __name__ == '__main__':
    main()

