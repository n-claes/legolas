class LegolasException(Exception):
    """
    Exception superclass to handle Legolas exceptions.

    Parameters
    ----------
    message : str
        The message to pass as error message.
    """

    def __init__(self, message=None):
        super().__init__(self, message)


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
