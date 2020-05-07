from argparse import ArgumentParser
from pathlib import Path
from post_processing import pylbo

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
    ps = pylbo.SingleSpectrum(ds, plot_continua=False)
    ps.plot_eigenfunctions()
    pylbo.show()

if __name__ == '__main__':
    _main()