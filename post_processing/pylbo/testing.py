import pickle


def pickle_dataseries_to_file(series, filepath):
    with open(filepath, "wb") as ostream:
        pickle.dump(series, ostream, pickle.HIGHEST_PROTOCOL)


def load_pickled_dataseries(filepath):
    with open(filepath, "rb") as istream:
        series = pickle.load(istream)
    return series
