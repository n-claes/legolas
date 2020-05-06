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
        return "There are still variables remaining in the supplied dictionary!\n" \
               "Check key names and values given. Remaining variables:\n" \
               "{}".format(self.file)


