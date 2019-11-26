import numpy as np
import os, sys

def check_file(filename):
    if not os.path.isfile(filename):
        print(">> file not found -> {}".format(filename))
        sys.exit(1)
    else:
        return


def read_stream_data(filename, content_type, rows, cols):
    """
    Method to read stream data from file. The parameters 'content_type', 'rows'
    and 'cols' must be supplied correctly in order to properly read the
    binary file with Numpy.
    :param filename: The name of the file (with path prepended if needed)
    :param content_type: The type of items in the file (double, complex, int...)
    :param rows: The number of rows in the file
    :param cols: The number of columns in the file
    """
    check_file(filename)

    nbitems = rows * cols
    print(">> Reading {}".format(filename))
    data = np.fromfile(filename, dtype=content_type, count=nbitems)
    results = data.reshape(rows, cols)

    return results

def read_config_file(filename):
    def read_bool(bool):
        return bool == 'T'

    check_file(filename)

    config_dict = {}

    configfile = open(filename, 'r')
    for line in configfile:
        line = line.split(':')
        line = [l.strip() for l in line]

        if len(line) == 1:
            continue

        if line[1] == 'T' or line[1] == 'F':
            config_dict[line[0]] = read_bool(line[1])
        else:
            try:
                if "." in line[1]:
                    config_dict[line[0]] = float(line[1])
                else:
                    config_dict[line[0]] = int(line[1])
            except ValueError:
                config_dict[line[0]] = line[1]

    return config_dict
