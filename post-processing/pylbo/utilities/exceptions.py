class LEGOLASException(Exception):
    def __init__(self, ds=None, message=None):
        super().__init__(self, message)
        self.ds = ds
