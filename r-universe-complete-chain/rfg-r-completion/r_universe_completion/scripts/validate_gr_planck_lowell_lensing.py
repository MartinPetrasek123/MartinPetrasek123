#!/usr/bin/env python3
"""Regression check for the executed GR Planck low-ell+lensing reference."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from gr_planck_lowell_lensing import calculate_loglikes


EXPECTED = {
    "planck_2018_lowl.TT": -1.1711737587329509e1,
    "planck_2018_lowl.EE": -1.9806344786661174e2,
    "planck_2018_lensing.clik": -4.3955688553898575,
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packages-path", type=Path, required=True, help="Cobaya external packages directory")
    arguments = parser.parse_args()
    actual = calculate_loglikes(arguments.packages_path)
    for component, expected in EXPECTED.items():
        value = actual[component]
        if not math.isclose(value, expected, rel_tol=5.0e-6, abs_tol=1.0e-10):
            raise AssertionError(f"{component}: got {value:.16e}, expected {expected:.16e}")
    print("Official Planck 2018 GR low-ell+lensing regression passed")


if __name__ == "__main__":
    main()
