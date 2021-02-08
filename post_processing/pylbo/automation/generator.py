import copy
import f90nml
from pathlib import Path
from pylbo.automation.defaults import namelist_items
from pylbo.utilities.toolbox import transform_to_list
from pylbo.utilities.logger import pylboLogger
from pylbo.exceptions import ParfileGenerationError


def _validate_basename(basename):
    """
    Validates the basename for a given parfile.

    Parameters
    ----------
    basename : str
        The basename for a parfile. If not given, defaults to "parfile".

    Returns
    -------
    basename : str
        The basename for a parfile.

    """
    if basename is None:
        basename = "parfile"
    return basename


def _validate_output_dir(output_dir):
    """
    Validates and returns the output directory for the parfiles.

    Parameters
    ----------
    output_dir : str, ~os.PathLike
        The output directory to store the parfiles in. If not given, defaults to
        the current working directory.

    Raises
    ------
    NotADirectoryError
        If the output directory is not found.

    Returns
    -------
    output : ~os.PathLike
        The resolved path to the output directory with "parfiles" appended.

    """
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
    """
    Handles parfile generation.

    Parameters
    ----------
    parfile_dict : dict
        Dictionary containing the keys to be placed in the parfile.
    basename : str
        The basename for the parfile, the `.par` suffix is added automatically and is
        not needed. If multiple parfiles are generated, these
        will be prepended by a 4-digit number (e.g. 0003myparfile.par).
        If not provided, the basename will default to `parfile`.
    output_dir : str, ~os.PathLike
        Output directory where the parfiles are saved, defaults to the current
        working directory if not specified. A subdirectory called `parfiles` will be
        created in which the parfiles will be saved.
    """

    def __init__(self, parfile_dict, basename=None, output_dir=None):
        self.parfile_dict = copy.deepcopy(parfile_dict)
        self.basename = _validate_basename(basename)
        self.output_dir = _validate_output_dir(output_dir)
        self.nb_runs = self.parfile_dict.pop("number_of_runs", 1)
        self.parfiles = []
        self.container = {}

    def _get_and_check_item(self, namelist, name, allowed_dtypes):
        """
        Does typechecking on the various dictionary keys supplied to the parfile
        generator. Pops the key from the dictionary.

        Parameters
        ----------
        namelist : str
            One of the namelists ("gridlist", "savelist", etc.)
        name : str
            The key to check.
        allowed_dtypes : class, tuple
            Allowed types for that particular key. Either a single value or a tuple.

        Raises
        ------
        TypeError
            If the types do not match, e.g. if "gridpoints" is specified as a float
            value when it should be an integer.

        Returns
        -------
        item : any
            The item popped from the dictionary corresponding to the given key.
        """
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
        """
        Creates one major namelist from the given dictionary.

        Raises
        ------
        ParfileGenerationError
            - If the original dictionary is not empty after everything should be popped
            - If there is an inconsistency between array sizes of dictionary items
        """
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
                # make sure that sigma has complex values
                if key == "sigma":
                    values_list = [complex(value) for value in values_list]
                self.container.get(namelist).update({key: values_list})

    def generate_parfiles(self):
        """
        Creates separate parfiles from the main namelist container and writes
        the individual parfiles to disk.

        Returns
        -------
        parfiles : list of str
            List containing the paths of the parfiles, can be passed to the legolas
            runner.

        """
        run_dict = {key: {} for key in self.container.keys()}
        # savelist must be present
        try:
            run_dict["savelist"]
        except KeyError:
            run_dict.update({"savelist": {}})

        for current_run in range(self.nb_runs):
            prefix = "{:04d}".format(current_run + 1)
            if self.nb_runs == 1:
                prefix = ""
            # generate dictionary for this specific run
            for namelist, items in self.container.items():
                for key, values in items.items():
                    run_dict[namelist].update({key: values[current_run]})
            # parfile name
            parfile_name = f"{prefix}{self.basename}.par"
            # datfile name (no extension .dat needed)
            datfile_name = "".join(
                [
                    f"{prefix}",
                    run_dict["savelist"].get("basename_datfile", self.basename),
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
