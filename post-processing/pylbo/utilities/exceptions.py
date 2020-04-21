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


