def custom_enumerate(iterable, start=0, step=1):
    for itr in iterable:
        yield start, itr
        start += step
