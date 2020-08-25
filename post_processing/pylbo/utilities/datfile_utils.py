import struct
import numpy as np

SIZE_CHAR = struct.calcsize('c')
SIZE_INT = struct.calcsize('i')
SIZE_BOOL = struct.calcsize('i')  # fortran logical is 4-byte integer
SIZE_DOUBLE = struct.calcsize('d')
SIZE_COMPLEX = struct.calcsize(2 * 'd')  # complex is 2 times double-byte

ALIGN = '='


def get_header(istream):
    """
    Retrieves the header from a given datfile.

    Parameters
    ----------
    istream : ~io.BufferedReader
        Datfile opened in binary mode.

    Returns
    -------
    h : dict
        Dictionary containing all header information.
    """
    istream.seek(0)
    h = {}

    # version information was added afterwards, check for this
    try:
        version_name = "datfile_version"
        fmt = ALIGN + len(version_name) * 'c'
        hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        assert b''.join(hdr).strip().decode() == version_name
        fmt = ALIGN + 'i'
        VERSION, = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    except AssertionError:
        istream.seek(0)
        VERSION = 0
    h['datfile_version'] = VERSION

    # read maximal string length and length of strings in arrays
    fmt = ALIGN + 2 * 'i'
    str_len, str_len_arr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))

    # read geometry
    fmt = ALIGN + str_len * 'c'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    h['geometry'] = b''.join(hdr).strip().decode()

    # read x_start, x_end and 4 different gridpoints (hence 2 doubles 'd' and 4 integers 'i')
    fmt = ALIGN + 2 * 'd' + 4 * 'i'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    h['x_start'], h['x_end'], h['gridpts'], h['gauss_gridpts'], h['matrix_gridpts'], h['ef_gridpts'] = hdr
    # read gamma
    fmt = ALIGN + 'd'
    h['gamma'], = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # read equilibrium type, this is a string of 'str_length'
    fmt = ALIGN + str_len * 'c'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    h['eq_type'] = b''.join(hdr).strip().decode()
    # read eigenfunctions boolean
    fmt = ALIGN + 'i'   # a fortran logical is a 4 byte integer
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))    # this is either (0,) or (1,) (F or T)
    h['eigenfuncs_written'] = bool(*hdr)    # bool casts 0 to False, everything else to True
    # # read matrices boolean
    fmt = ALIGN + 'i'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    h['matrices_written'] = bool(*hdr)

    # read number of parameters saved
    fmt = ALIGN + 'i'
    nb_params, = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # read parameter names, there are 'nb_params' of them each with size str_len_arr
    fmt = ALIGN + nb_params * str_len_arr * 'c'
    param_names = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    param_names = [b''.join(param_names[i:i + str_len_arr]).strip().decode()
                   for i in range(0, len(param_names), str_len_arr)]
    # read parameter values, these are all floats (including NaN)
    fmt = ALIGN + nb_params * 'd'
    param_values = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    params = {}
    for name, value in zip(param_names, param_values):
        # skip NaN parameter values, these are not used
        if np.isnan(value):
            continue
        params[name] = value
    h['params'] = params

    # read names for saved equilibrium arrays
    fmt = ALIGN + 'i'
    nb_equilvars, = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # read equilibrium array names, 'nb_equilvars' each with size str_len_arr (see mod_output)
    fmt = ALIGN + nb_equilvars * str_len_arr * 'c'
    equil_names = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    equil_names = [b''.join(equil_names[i:i+str_len_arr]).strip().decode()
                   for i in range(0, len(equil_names), str_len_arr)]
    h['equil_names'] = equil_names

    # units-related stuff
    units = {}
    fmt = ALIGN + 'i'
    h['cgs'] = bool(struct.unpack(fmt, istream.read(struct.calcsize(fmt))))
    fmt = ALIGN + 11 * 'd'  # there are 11 unit normalisations in total
    unit_values = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    unit_names = ['unit_length', 'unit_time', 'unit_density', 'unit_velocity', 'unit_temperature', 'unit_pressure',
                  'unit_magneticfield', 'unit_numberdensity', 'unit_luminosity', 'unit_conduction', 'unit_resistivity']
    for name, value in zip(unit_names, unit_values):
        units[name] = value
    h['units'] = units

    # @Devnote Niels: do not read in data here but simply determine the data offset in the binary file.
    # These offsets are later used to actually read in the data and immediately set it as class attributes.
    # This prevents double loading.

    # eigenvalue offset
    byte_size = h['matrix_gridpts'] * SIZE_COMPLEX  # eigenvalues are complex
    offsets = {'eigenvalues': istream.tell()}
    istream.seek(istream.tell() + byte_size)
    # grid offset
    byte_size = h['gridpts'] * SIZE_DOUBLE
    offsets.update({'grid': istream.tell()})
    istream.seek(istream.tell() + byte_size)
    # grid_gauss offset
    byte_size = h['gauss_gridpts'] * SIZE_DOUBLE
    offsets.update({'grid_gauss': istream.tell()})
    istream.seek(istream.tell() + byte_size)
    # equilibrium arrays offset
    byte_size = h['gauss_gridpts'] * len(equil_names) * SIZE_DOUBLE
    offsets.update({'equil_arrays': istream.tell()})
    istream.seek(istream.tell() + byte_size)

    # if eigenfunctions are written, read names and include offsets
    if h['eigenfuncs_written']:
        # read eigenfunction names
        fmt = ALIGN + 'i'
        nb_eigenfuncs, = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        fmt = ALIGN + nb_eigenfuncs * str_len_arr * 'c'
        ef_names = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        ef_names = [b''.join(ef_names[i:i + str_len_arr]).strip().decode()
                    for i in range(0, len(ef_names), str_len_arr)]
        h['ef_names'] = ef_names
        # ef_grid offset
        byte_size = h['ef_gridpts'] * SIZE_DOUBLE
        offsets.update({'ef_grid': istream.tell()})
        istream.seek(istream.tell() + byte_size)
        # eigenfunction offsets
        byte_size = h['ef_gridpts'] * h['matrix_gridpts'] * nb_eigenfuncs * SIZE_COMPLEX
        offsets.update({'ef_arrays': istream.tell()})
        istream.seek(istream.tell() + byte_size)

    # if matrices are written, include amount of nonzero elements and offsets
    if h['matrices_written']:
        fmt = ALIGN + 2 * 'i'
        hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        h['nonzero_B_elements'], h['nonzero_A_elements'] = hdr

        # matrix B offset: this is written as (row, column, value)
        byte_size = (2 * SIZE_INT + SIZE_DOUBLE) * h['nonzero_B_elements']
        offsets.update({'matrix_B': istream.tell()})
        istream.seek(istream.tell() + byte_size)
        # matrix A offset
        offsets.update({'matrix_A': istream.tell()})
    h['offsets'] = offsets
    return h


def read_grid(istream, header):
    """
    Retrieves the base grid from the datfile.

    Parameters
    ----------
    istream : ~io.BufferedReader
        Datfile opened in binary mode.
    header : dict
        Dictionary containing the datfile header, output from :func:`get_header`.

    Returns
    -------
    grid : numpy.ndarray(dtype=float, ndim=1)
        The base grid from the datfile.
    """
    istream.seek(header['offsets']['grid'])
    fmt = ALIGN + header['gridpts'] * 'd'
    grid = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    return np.asarray(grid)


def read_grid_gauss(istream, header):
    """
    Retrieves the Gaussian grid from the datfile.

    Parameters
    ----------
    istream : ~io.BufferedReader
        Datfile opened in binary mode.
    header : dict
        Dictionary containing the datfile header, output from :func:`get_header`.

    Returns
    -------
    grid_gauss : numpy.ndarray(dtype=float, ndim=1)
        The Gaussian grid from the datfile.
    """
    istream.seek(header['offsets']['grid_gauss'])
    fmt = ALIGN + header['gauss_gridpts'] * 'd'
    grid_gauss = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    return np.asarray(grid_gauss)


def read_ef_grid(istream, header):
    """
    Retrieves the eigenfunction grid from the datfile.

    Parameters
    ----------
    istream : ~io.BufferedReader
        Datfile opened in binary mode.
    header : dict
        Dictionary containing the datfile header, output from :func:`get_header`.

    Returns
    -------
    ef_grid : numpy.ndarray(dtype=float, ndim=1)
        The eigenfunction grid from the datfile.
    """
    istream.seek(header['offsets']['ef_grid'])
    fmt = ALIGN + header['ef_gridpts'] * 'd'
    ef_grid = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    return np.asarray(ef_grid)


def read_eigenvalues(istream, header, omit_large_evs=True):
    """
    Reads the eigenvalues from the datfile, and optionally omits
    the large (larger than `1e15`) values.

    Parameters
    ----------
    istream : ~io.BufferedReader
        Datfile opened in binary mode.
    header : dict
        Dictionary containing the datfile header, output from :func:`get_header`.
    omit_large_evs : bool
        If `True`, all eigenvalues with a modulus larger than `1e15` are omitted.

    Returns
    -------
    eigenvalues : numpy.ndarray(dtype=complex, ndim=1)
        The eigenvalues from the datfile, with optionally omitted large values.
    """
    istream.seek(header['offsets']['eigenvalues'])
    fmt = ALIGN + 2 * header['matrix_gridpts'] * 'd'
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    # hdr is a 1D list with [real, imag, real, imag, real, imag...]
    reals = hdr[::2]
    imags = hdr[1::2]
    eigenvalues = np.asarray([complex(x, y) for x, y in zip(reals, imags)])
    if omit_large_evs:
        eigenvalues = eigenvalues[np.where(np.absolute(eigenvalues) < 1e15)]
    return eigenvalues


def read_equilibrium_arrays(istream, header):
    """
    Reads the equilibrium arrays from the datfile.

    Parameters
    ----------
    istream : ~io.BufferedReader
        Datfile opened in binary mode.
    header : dict
        Dictionary containing the datfile header, output from :func:`get_header`.

    Returns
    -------
    equil_arrays : dict
        Dictionary containing the equilibrium arrays, with keys given by
        `header['equil_names']`.
    """
    istream.seek(header['offsets']['equil_arrays'])
    equil_arrays = {}
    for name in header['equil_names']:
        fmt = ALIGN + header['gauss_gridpts'] * 'd'
        equil_array = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        equil_arrays.update({name: np.asarray(equil_array)})
    return equil_arrays


def read_eigenfunction(istream, header, ef_index):
    """
    Reads a single eigenfunction from the datfile. Eigenfunctions are read in on-the-fly,
    to prevent having to load the entire array into memory (which can quickly be a few
    100 Mb for larger datasets).

    Parameters
    ----------
    istream : ~io.BufferedReader
        Datfile opened in binary mode.
    header : dict
        Dictionary containing the datfile header, output from :func:`get_header`.
    ef_index : int
        The index of the eigenfunction in the matrix. This value corresponds to the
        index of the accompanying eigenvalue in the
        :attr:`~pylbo.LegolasDataContainer.eigenvalues` attribute.
        These indices can be retrieved through
        :func:`~pylbo.LegolasDataContainer.get_nearest_eigenvalues`.

    Returns
    -------
    eigenfunctions : dict
        A dictionary containing the eigenfunctions for all variables with keys given
        by the names of `header['ef_names']`. The eigenfunctions correspond to a specific
        eigenvalue, associated with the same `ef_index`.
    """
    ef_offset = header['offsets']['ef_arrays']
    ef_gridpts = header['ef_gridpts']
    matrix_gridpts = header['matrix_gridpts']
    eigenfunctions = {}

    # Fortran writes in column-major order, meaning column per column.
    # This makes it quite convenient to extract a single eigenfunction from the datfile
    matrix_bytesize = ef_gridpts * matrix_gridpts * SIZE_COMPLEX
    ef_bytesize = ef_gridpts * SIZE_COMPLEX
    for name in header['ef_names']:
        name_idx = header['ef_names'].index(name)
        # move pointer to correct place of current name matrix in datfile
        istream.seek(ef_offset + name_idx * matrix_bytesize)
        # move pointer to requested eigenfunction
        # 'ef_index' here is the index of the corresponding eigenvalue in its array
        istream.seek(istream.tell() + ef_index * ef_bytesize)
        # read in single eigenfunction
        fmt = ALIGN + ef_gridpts * 2 * 'd'
        hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        reals = hdr[::2]
        imags = hdr[1::2]
        ef_values = np.asarray([complex(x, y) for x, y in zip(reals, imags)])
        eigenfunctions.update({name: ef_values})
    return eigenfunctions


def read_matrix_B(istream, header):
    """
    Reads the B-matrix from the datfile.

    Parameters
    ----------
    istream : ~io.BufferedReader
        Datfile opened in binary mode.
    header : dict
        Dictionary containing the datfile header, output from :func:`get_header`.

    Returns
    -------
    rows : numpy.ndarray(dtype=int, ndim=1)
        Array containing the row indices of the non-zero matrix entries.
    cols : numpy.ndarray(dtype=int, ndim=1)
        Array containing the column indices of the non-zero matrix entries.
    vals : numpy.ndarray(dtype=float, ndim=1)
        Array containing the non-zero B-matrix values corresponding to
        the rows and column indices.
    """
    istream.seek(header['offsets']['matrix_B'])
    fmt = ALIGN + (2 * 'i' + 'd') * header['nonzero_B_elements']
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    rows = np.asarray(hdr[::3])     # rows are 1, 4, 7, 10 etc (Fortran indexing)
    cols = np.asarray(hdr[1::3])    # columns are 2, 5, 8, 11 etc (Fortran indexing)
    vals = np.asarray(hdr[2::3])    # values are 3, 6, 9, 12 etc (Fortran indexing)
    return rows, cols, vals


def read_matrix_A(istream, header):
    """
    Reads the A-matrix from the datfile.

    Parameters
    ----------
    istream : ~io.BufferedReader
        Datfile opened in binary mode.
    header : dict
        Dictionary containing the datfile header, output from :func:`get_header`.

    Returns
    -------
    rows : numpy.ndarray(dtype=int, ndim=1)
        Array containing the row indices of the non-zero matrix entries.
    cols : numpy.ndarray(dtype=int, ndim=1)
        Array containing the column indices of the non-zero matrix entries.
    vals : numpy.ndarray(dtype=complex, ndim=1)
        Array containing the non-zero A-matrix values corresponding to
        the row and column indices.
    """
    istream.seek(header['offsets']['matrix_A'])
    fmt = ALIGN + (2 * 'i' + 2 * 'd') * header['nonzero_A_elements']
    hdr = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
    rows = np.asarray(hdr[::4])
    cols = np.asarray(hdr[1::4])
    reals = hdr[2::4]
    imags = hdr[3::4]
    vals = np.asarray([complex(x, y) for x, y in zip(reals, imags)])
    return rows, cols, vals
