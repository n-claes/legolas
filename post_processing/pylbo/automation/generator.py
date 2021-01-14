import numpy as np
import f90nml
from pathlib import Path
from pylbo.automation.defaults import namelist_items
from pylbo.utilities.toolbox import transform_to_list
from pylbo.utilities.logger import pylboLogger
from pylbo.exceptions import ParfileGenerationError


def _validate_basename(basename):
    if basename is None:
        basename = "parfile"
    return basename


def _validate_output_dir(output_dir):
    if output_dir is None:
        output_dir = Path.cwd()
    output = Path(output_dir).resolve()
    if not output.is_dir():
        raise NotADirectoryError(output)
    output = (output / "parfiles").resolve()
    if not output.is_dir():
        Path.mkdir(output)
    return output


class ParfileGenerator:
    def __init__(self, parfile_dict, basename=None, output_dir=None):
        self.parfile_dict = parfile_dict
        self.basename = _validate_basename(basename)
        self.output_dir = _validate_output_dir(output_dir)
        self.nb_runs = self.parfile_dict.pop("number_of_runs", 1)
        self.parfiles = []
        self.container = {}

    def _get_and_check_item(self, namelist, name, allowed_dtypes):
        item = self.parfile_dict.pop(name, None)
        if item is not None:
            item_aslist = transform_to_list(item)
            for obj in item_aslist:
                if not isinstance(obj, allowed_dtypes):
                    raise TypeError(
                        f"namelist '{namelist}' expected item '{name}' to be of "
                        f"type {allowed_dtypes} but got {type(obj)}. \n"
                        f"item value = {item}"
                    )
        return item

    def create_namelist_from_dict(self):
        for namelist, items in namelist_items.items():
            # update container
            self.container.update({namelist: {}})
            # loop over names for this specific namelist
            for name, dtypes in items:
                obj = self._get_and_check_item(namelist, name, dtypes)
                if obj is not None:
                    obj = transform_to_list(obj)
                    self.container[namelist].update({name: obj})
            # if this namelist is still empty, remove it
            if self.container[namelist] == {}:
                self.container.pop(namelist)
        # account for parameters separately
        params = self.parfile_dict.pop("parameters", {})
        if len(params) != 0:
            self.container["equilibriumlist"].update({"use_defaults": [False]})
            for name, param in params.items():
                params.update({name: transform_to_list(param)})
        self.container.update({"paramlist": params})

        # here the original dictionary should be empty,
        # something went wrong if it isn't
        if len(self.parfile_dict) != 0:
            raise ParfileGenerationError(self.parfile_dict)

        # update container for number of runs, all items are lists
        for namelist, items in self.container.items():
            for key, values in items.items():
                if len(values) == 1:
                    values_list = values * self.nb_runs
                elif len(values) == self.nb_runs:
                    values_list = values
                else:
                    raise ParfileGenerationError(items, self.nb_runs, key)
                self.container.get(namelist).update({key: values_list})

    def generate_parfiles(self):
        run_dict = {key: {} for key in self.container.keys()}
        # savelist must be present
        try:
            run_dict["savelist"]
        except KeyError:
            run_dict.update({"savelist", {}})

        for current_run in range(self.nb_runs):
            prefix = "{:04d}".format(current_run + 1)
            if self.nb_runs == 1:
                prefix = ""
            # generate dictionary for this specific run
            for namelist, items in self.container.items():
                for key, values in items.items():
                    run_dict[namelist].update({key: values[current_run - 1]})
            # parfile name
            parfile_name = f"{prefix}{self.basename}.par"
            # datfile name (no extension .dat needed)
            datfile_name = "".join(
                [
                    f"{prefix}",
                    run_dict["savelist"].get("basename_datfile", self.basename)
                ]
            )
            run_dict["savelist"].update({"basename_datfile": datfile_name})
            # logfile name (no extension .log needed)
            logfile_name = run_dict["savelist"].get("basename_logfile", None)
            if logfile_name is not None:
                logfile_name = f"{prefix}{logfile_name}"
                run_dict["savelist"].update({"basename_logfile": logfile_name})

            # set paths and write parfile
            parfile_path = (self.output_dir / parfile_name).resolve()
            self.parfiles.append(str(parfile_path))
            f90nml.write(run_dict, parfile_path, force=True)
        pylboLogger.info(f"parfiles generated and saved to {self.output_dir}")
        return self.parfiles
