from argparse import ArgumentParser
from pathlib import Path
import pylbo
from pylbo.utilities.exceptions import EigenfunctionsNotPresent, MatricesNotPresent


def _main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--datfile', dest='datfile')
    args = parser.parse_args()
    datfile = args.datfile

    if datfile is None:
        datfile = Path('output/datfile.dat').resolve()
        if not datfile.is_file():
            raise FileNotFoundError(datfile)

    ds = pylbo.load(datfile)
    ps = pylbo.SingleSpectrum(ds)
    ps.plot_spectrum(annotate_continua=True)
    try:
        ps.plot_eigenfunctions(merge_figs=True)
    except EigenfunctionsNotPresent:
        pass
    if ds.gridpts < 20:
        try:
            pylbo.plot_matrices(ds)
        except MatricesNotPresent:
            pass
    ps.plot_equilibria()
    pylbo.show()


if __name__ == '__main__':
    _main()
