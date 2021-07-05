import copy
import shutil
from pathlib import Path
from argparse import ArgumentParser
import pylbo
import pylbo.testing

from regression_tests.regression import tests_to_run, multirun_tests


def overwrite_files(datfile):
    if Path.is_file(datfile):
        print(f"{datfile.name} already exists!")
        return input("overwrite? ").lower() in ("yes", "y")
    return True


def single_files():
    parfiles = []
    names = []
    for test in tests_to_run:
        name = test["name"]
        print("=" * 50)
        print(">> generating {}".format(name))
        config = copy.deepcopy(test["config"])
        output_folder = str(test["answer_datfile"].parent)
        config.update(
            {
                "basename_datfile": test["answer_datfile"].stem,
                "basename_logfile": test["answer_logfile"].stem,
                "output_folder": output_folder,
            }
        )
        if overwrite_files(test["answer_datfile"]):
            parfile = pylbo.generate_parfiles(
                parfile_dict=config,
                basename=test["answer_datfile"].stem,
                output_dir=output_folder,
            )
            parfiles.append(*parfile)
            names.append(test["answer_datfile"].stem)
        else:
            print(f"Skipping {name}")
    if parfiles:
        print("=" * 50)
        print(
            f"Generator will create or overwrite the "
            f"following files: {[name for name in names]}"
        )
        force = input("Are you sure? ")
        if force.lower() in ("yes", "y"):
            pylbo.run_legolas(parfiles, remove_parfiles=True, nb_cpus=4)


def multirun_files():
    for test in multirun_tests:
        name = test["name"]
        print("=" * 50)
        print(">> generating {}".format(name))
        config = copy.deepcopy(test["config"])
        pickled_file = test["answer_series"]
        if overwrite_files(pickled_file):
            # create temporary folder
            tmp = (test["answer_series"].parent / "tmp").resolve()
            if tmp.is_dir():
                shutil.rmtree(tmp)
            tmp.mkdir()
            config["output_folder"] = str(tmp)
            config["basename_datfile"] = f"{name}_tmp"
            parfiles = pylbo.generate_parfiles(
                parfile_dict=config, basename=pickled_file.stem, output_dir=tmp
            )
            print(
                f"Generator will create or overwrite the "
                f"following file: {pickled_file}"
            )
            if input("Are you sure? ").lower() in ("yes", "y"):
                pylbo.run_legolas(parfiles, remove_parfiles=True, nb_cpus=4)
                # load file, pickle, save, remove temporary folder
                series = pylbo.load_series(sorted(tmp.glob(f"*{name}_tmp.dat")))
                pylbo.testing.pickle_dataseries_to_file(series, pickled_file)
            shutil.rmtree(tmp)
        else:
            print(f"Skipping {name}")


def main():
    parser = ArgumentParser()
    parser.add_argument("--multi", action="store_true")
    args = parser.parse_args()

    if args.multi:
        multirun_files()
    else:
        single_files()


if __name__ == "__main__":
    main()
