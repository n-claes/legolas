class LegolasException(Exception):
    """
    Exception superclass to handle Legolas exceptions.

    Parameters
    ----------
    ds : ~pylbo.LegolasDataContainer
        The :class:`~pylbo.LegolasDataContainer` instance currently loaded.
    message : str
        The message to pass as error message.
    """

    def __init__(self, ds=None, message=None):
        super().__init__(self, message)
        self.ds = ds


class InvalidLegolasFile(LegolasException):
    """
    Handles trying to load invalid Legolas files.

    Parameters
    ----------
    file : str or ~os.PathLike
        The path to the file.
    """

    def __init__(self, file):
        self.file = file

    def __str__(self):
        return f"{self.file} is not a valid Legolas file"


class EigenfunctionsNotPresent(LegolasException):
    """
    Handles trying to query for eigenfunctions when these
    are not present in the datfile.

    Parameters
    ----------
    file : str or ~os.PathLike
        The path to the file.
    """

    def __init__(self, file):
        self.file = file

    def __str__(self):
        return f"No eigenfunctions present in {self.file}"


class MatricesNotPresent(LegolasException):
    """
    Handles trying to query for matrices when these
    are not present in the datfile.

    Parameters
    ----------
    file : str or ~os.PathLike
        The path to the file.
    """

    def __init__(self, file):
        self.file = file

    def __str__(self):
        return f"No matrices present in {self.file}"


class InconsistentMultirunFile(LegolasException):
    """
    Throws an error when loading multiple datfiles, and one
    or more contain a different equilibrium type.

    Parameters
    ----------
    file : str or ~os.PathLike
        The path to the file.
    found : str
        The name of the inconsistent equilibrium type encountered.
    expected : str
        The name of the expected equilibrium type.
    """

    def __init__(self, file, found, expected):
        self.file = file
        self.found = found
        self.expected = expected

    def __str__(self):
        return (
            f"Different equilibrium encountered in {self.file}. \n"
            f"Expected: {self.found}\n"
            f"Found   : {self.expected}"
        )


class DictNotEmpty(LegolasException):
    """
    Throws an error during parfile generation, when the supplied dictionary
    has not correctly been popped. This can occur due to a typing error in the keys,
    or supplying a key that can not be supplied.

    Parameters
    ----------
    file : str or ~os.PathLike
        The path to the file.
    """

    def __init__(self, file):
        self.file = file

    def __str__(self):
        return (
            f"There are still variables remaining in the supplied dictionary "
            f"after parfile generation!\n"
            f"Check key names and values given. Remaining variables:\n"
            f"{self.file}"
        )


class UnknownPrecodedRun(LegolasException):
    """
    Throws an error when trying to obtain a precoded dictionary
    with an unknown name.

    Parameters
    ----------
    given_name : str
        The name of the key to be obtained.
    available_names : list
        A list of available keys to choose from.
    """

    def __init__(self, given_name, available_names):
        self.given_name = given_name
        self.available_names = list(available_names)

    def __str__(self):
        return (
            f"Unknown precoded run: {self.given_name}\n"
            "Choose from {self.available_names}"
        )


class InconsistentNumberOfRuns(LegolasException):
    """
    Throws an error during parfile generation, when there is
    an inconsistency between the number of runs given and the size of
    one or more arrays in the configuration dictionary.

    Parameters
    ----------
    nb_runs : int
        The number of runs given.
    key : str
        The key of the item for which the inconsistency occurred
    par_dict : dict
        The configuration dictionary supplied.
    """

    def __init__(self, nb_runs, key, par_dict):
        self.nb_runs = nb_runs
        self.key = key
        self.par_dict = par_dict

    def __str__(self):
        return (
            f"Key 'number_of_runs' is inconsistent with at least one "
            f"supplied parameter.\nNumber of runs: {self.nb_runs}\n"
            f"Length of {self.key}: {self.par_dict[self.key]}"
        )
