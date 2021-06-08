import struct
import numpy as np
import matplotlib.pyplot as plt
from pylbo.utilities.datfile_utils import ALIGN


def get_debug_coolingcurves(filename):
    coolcurves = {}
    with open(filename, "rb") as istream:
        istream.seek(0)
        fmt = ALIGN + "i"
        (size_tables,) = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        fmt = ALIGN + size_tables * "d"
        coolcurves["T_table"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        coolcurves["L_table"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        fmt = ALIGN + "i"
        (size_curves,) = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        fmt = ALIGN + size_curves * "d"
        coolcurves["T_curve"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        coolcurves["L_curve"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        coolcurves["dLdT_curve"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
    return coolcurves


def get_debug_atmocurves(filename):
    atmocurves = {}
    with open(filename, "rb") as istream:
        istream.seek(0)
        fmt = ALIGN + "i"
        (size_tables,) = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        fmt = ALIGN + size_tables * "d"
        atmocurves["h_table"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        atmocurves["T_table"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        atmocurves["nh_table"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        fmt = ALIGN + "i"
        (size_curves,) = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
        fmt = ALIGN + size_curves * "d"
        atmocurves["h_curve"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        atmocurves["T_curve"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        atmocurves["nh_curve"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
        atmocurves["dTdr_curve"] = np.asarray(
            struct.unpack(fmt, istream.read(struct.calcsize(fmt))),
            dtype=float,
        )
    return atmocurves


def plot_debug_coolingcurves(filename):
    coolcurves = get_debug_coolingcurves(filename)

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(10, 6))
    axl.loglog(coolcurves["T_table"], coolcurves["L_table"], ".b", label="table")
    axl.loglog(coolcurves["T_curve"], coolcurves["L_curve"], "-r", label="curve")
    axl.set_xlabel("temperature [K]")
    axl.set_ylabel(r"$\Lambda(T)$ [erg s$^{-1}$ cm$^3$]")
    axl.legend(loc="best")
    axr.semilogx(
        coolcurves["T_curve"],
        coolcurves["dLdT_curve"],
        color="red",
        label="legolas 6th O",
    )
    axr.semilogx(
        coolcurves["T_curve"],
        np.gradient(coolcurves["L_curve"], coolcurves["T_curve"], edge_order=2),
        label="numpy 2nd O",
    )
    axr.set_xlabel("temperature [K]")
    axr.set_ylabel(r"$\frac{d\Lambda(T)}{dT}$")
    axr.legend(loc="best")
    fig.suptitle(f"cooling curves, interpolated at {len(coolcurves['L_curve'])} points")
    fig.tight_layout()
    return (fig, (axl, axr))


def plot_debug_atmocurves(filename):
    atmocurves = get_debug_atmocurves(filename)

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(10, 6))
    axl.semilogy(atmocurves["h_table"], atmocurves["T_table"], ".b", label="table")
    axl.semilogy(atmocurves["h_curve"], atmocurves["T_curve"], "-r", label="curve")
    axl.set_xlabel("height [cm]")
    axl.set_ylabel("temperature [K]")
    axl.legend(loc="best")
    axr.plot(
        atmocurves["h_curve"],
        atmocurves["dTdr_curve"],
        color="red",
        label="legolas 6th O",
    )
    axr.plot(
        atmocurves["h_curve"],
        np.gradient(atmocurves["T_curve"], atmocurves["h_curve"], edge_order=2),
        label="numpy 2nd O",
    )
    axr.set_xlabel("height [cm]")
    axr.set_ylabel(r"$\frac{dT}{dr}$")
    axr.legend(loc="best")
    fig.suptitle(
        f"atmosphere curves, interpolated at {len(atmocurves['T_curve'])} points"
    )
    fig.tight_layout()
    return (fig, (axl, axr))
