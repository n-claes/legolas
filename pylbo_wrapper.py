from argparse import ArgumentParser
from pathlib import Path
import pylbo


def _main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--datfile", dest="datfile")
    args = parser.parse_args()
    datfile = args.datfile

    if datfile is None:
        datfile = Path("output/datfile.dat").resolve()
        if not datfile.is_file():
            raise FileNotFoundError(datfile)

    ds = pylbo.load(datfile)
    pylbo.plot_equilibrium(ds)
    p = pylbo.plot_spectrum(ds)
    p.add_continua()
    if ds.header["eigenfuncs_written"]:
        p.add_eigenfunctions()
    p.show()


if __name__ == "__main__":
    _main()
