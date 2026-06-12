#!/usr/bin/env python3

import os
import re
from pathlib import Path

import numpy as np


LEGOLAS_FILE = (
    Path(os.environ["LEGOLASDIR"])
    / "src/physics/cooling_curves/mod_radloss_tables.f08"
)

AMRVAC_FILE = (
    Path(os.environ["AMRVAC_DIR"])
    / "src/physics/mod_radloss_tables.t"
)


def read_file(fname):
    text = Path(fname).read_text()

    # Remove comments
    text = re.sub(r"!.*$", "", text, flags=re.MULTILINE)

    # Remove continuation markers
    text = text.replace("&", " ")

    return text


def parse_data_blocks(text):
    """
    Parse Fortran DATA statements.

    Returns
    -------
    dict[str, np.ndarray]
    """

    pattern = re.compile(
        r"data\s+([A-Za-z0-9_]+)\s*/(.*?)/",
        re.IGNORECASE | re.DOTALL,
    )

    number_pattern = re.compile(
        r"""
        [-+]?
        (?:
            \d+\.\d* |
            \.\d+    |
            \d+
        )
        (?:[DdEe][-+]?\d+)?
        (?:_[A-Za-z0-9_]+)?
        """,
        re.VERBOSE,
    )

    tables = {}

    for name, values in pattern.findall(text):

        numbers = []

        for token in number_pattern.findall(values):

            token = re.sub(r"_[A-Za-z0-9_]+$", "", token)
            token = token.replace("D", "E").replace("d", "e")

            try:
                numbers.append(float(token))
            except ValueError:
                raise ValueError(
                    f"Failed to parse token '{token}' "
                    f"in table '{name}'"
                )

        tables[name] = np.asarray(numbers)

    return tables


def translate_amrvac_names(amrvac_tables):
    """
    Convert AMRVAC names to LEGOLAS names.
    """

    translated = {}

    for name, arr in amrvac_tables.items():

        if name.startswith("t_"):
            translated["logT_" + name[2:]] = arr

        elif name.startswith("l_"):
            translated["logL_" + name[2:]] = arr

        else:
            translated[name] = arr

    return translated


def compare_tables(legolas, amrvac, rtol=1e-12, atol=1e-12):

    common = sorted(set(legolas) & set(amrvac))

    print(f"\nFound {len(common)} common tables\n")

    n_ok = 0
    n_bad = 0

    for name in common:

        a = legolas[name]
        b = amrvac[name]

        if len(a) != len(b):

            n_bad += 1

            print(
                f"{name}: LENGTH MISMATCH "
                f"(LEGOLAS={len(a)}, AMRVAC={len(b)})"
            )

            continue

        if np.allclose(a, b, rtol=rtol, atol=atol):

            n_ok += 1

        else:

            n_bad += 1

            diff = np.abs(a - b)
            idx = np.argmax(diff)

            print(f"\n{name}: VALUES DIFFER")

            print(f"  index      : {idx}")
            print(f"  LEGOLAS    : {a[idx]:.16e}")
            print(f"  AMRVAC     : {b[idx]:.16e}")
            print(f"  abs diff   : {diff[idx]:.16e}")

    print("\nSummary")
    print("-------")
    print(f"Matching tables      : {n_ok}")
    print(f"Problematic tables   : {n_bad}")

    missing_legolas = sorted(set(amrvac) - set(legolas))
    missing_amrvac = sorted(set(legolas) - set(amrvac))

    if missing_legolas:
        print("\nMissing in LEGOLAS:")
        for name in missing_legolas:
            print("   ", name)

    if missing_amrvac:
        print("\nMissing in AMRVAC:")
        for name in missing_amrvac:
            print("   ", name)


def main():

    print(f"Reading {LEGOLAS_FILE}")
    print(f"Reading {AMRVAC_FILE}")

    legolas = parse_data_blocks(read_file(LEGOLAS_FILE))
    amrvac = parse_data_blocks(read_file(AMRVAC_FILE))

    amrvac = translate_amrvac_names(amrvac)

    compare_tables(legolas, amrvac)


if __name__ == "__main__":
    main()