import pylbo
import copy
from pathlib import Path
from regression_tests.regression import tests_to_run


def overwrite_files(datfile):
    if Path.is_file(datfile):
        print(f"{datfile.name} already exists!")
        force = input("overwrite? ")
        if force.lower() in ("yes", "y"):
            return True
        else:
            return False
    return True


def main():
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
        else:
            pass


if __name__ == "__main__":
    main()
