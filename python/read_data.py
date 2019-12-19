import numpy as np
import os, sys
try:
    import f90nml
except ImportError:
    print('f90nml package has to be installed first. \n'
          '>> conda install -c conda-forge f90nml')
    sys.exit(1)


def check_file(filename):
    """
    Checks if the provided file exists
    :param filename: path to file
    :return: true if path exists, false otherwise
    """
    if not os.path.isfile(filename):
        print(">> file not found -> {}".format(filename))
        sys.exit(1)


def read_config_file(config_file):
    """
    Reads the configuration namelist and processes it as a nested dictionary
    :param config_file: namelist to read
    :return: nested dictionary, access through eg. config_file['gridlist']['geometry']
    """
    check_file(config_file)

    config_dict = f90nml.read(config_file)

    # check strings for leading or trailing spaces and remove them
    # iterate over different namelists
    for namelist in config_dict:
        # iterate over keys in each namelist
        for namelist_key in config_dict[namelist]:
            # if corresponding value is a string, remove spaces
            if isinstance(config_dict[namelist][namelist_key], str):
                config_dict[namelist][namelist_key] = config_dict[namelist][namelist_key].strip()

    return config_dict


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
    stream_data = data.reshape(rows, cols)

    return stream_data
