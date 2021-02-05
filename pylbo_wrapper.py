from argparse import ArgumentParser
from pathlib import Path
import pylbo


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
    p2 = pylbo.plot_equilibrium(ds)
    p = pylbo.plot_spectrum(ds)
    p.add_continua()
    p.add_eigenfunctions()
    p.showall()


if __name__ == '__main__':
    _main()
