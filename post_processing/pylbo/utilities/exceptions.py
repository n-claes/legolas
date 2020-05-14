class LegolasException(Exception):
    def __init__(self, ds=None, message=None):
        super().__init__(self, message)
        self.ds = ds


class InvalidLegolasFile(LegolasException):
    def __init__(self, file):
        self.file = file

    def __str__(self):
        return "{} is not a valid Legolas file".format(self.file)


class EigenfunctionsNotPresent(LegolasException):
    def __init__(self, file):
        self.file = file

    def __str__(self):
        return "No eigenfunctions present in {}".format(self.file)


class MatricesNotPresent(LegolasException):
    def __init__(self, file):
        self.file = file

    def __str__(self):
        return "No matrices present in {}".format(self.file)


class InconsistentMultirunFile(LegolasException):
    def __init__(self, file, found, expected):
        self.file = file
        self.found = found
        self.expected = expected

    def __str__(self):
        return "Different equilibrium encountered in {}. \n" \
               "Expected: {}\n" \
               "Found   : {}".format(self.file, self.found, self.expected)


class DictNotEmpty(LegolasException):
    def __init__(self, file):
        self.file = file

    def __str__(self):
        return "There are still variables remaining in the supplied dictionary after parfile generation!\n" \
               "Check key names and values given. Remaining variables:\n" \
               "{}".format(self.file)


class UnknownPrecodedRun(LegolasException):
    def __init__(self, given_name, available_names):
        self.given_name = given_name
        self.available_names = list(available_names)

    def __str__(self):
        return "Unknown precoded run: {}\n" \
               "Choose from {}".format(self.given_name, self.available_names)


class InconsistentNumberOfRuns(LegolasException):
    def __init__(self, nb_runs, key, par_dict):
        self.nb_runs = nb_runs
        self.key = key
        self.par_dict = par_dict

    def __str__(self):
        return "Key 'number_of_runs' is inconsistent with at least one supplied parameter.\n" \
               "Number of runs: {}\n" \
               "Length of {}: {}".format(self.nb_runs, self.key, len(self.par_dict[self.key]))
