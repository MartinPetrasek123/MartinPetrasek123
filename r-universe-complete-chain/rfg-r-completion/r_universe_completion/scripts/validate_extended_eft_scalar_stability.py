#!/usr/bin/env python3
"""Regression test for the RFG-R extended-EFT scalar degeneracy result."""

from __future__ import annotations

import math

import numpy as np

from extended_eft_scalar_stability import scalar_stability_row
from rfg_regularized import RFGRegularizedParams


def main() -> None:
    params = RFGRegularizedParams()
    params.validate()
    rows = [
        scalar_stability_row(float(a), float(k), params)
        for a in np.logspace(-7.0, 0.0, 49)
        for k in np.logspace(-4.0, 5.0, 49)
    ]
    assert len(rows) == 49 * 49
    assert all(row["scalar_kinetic_degenerate"] == 1 for row in rows)
    assert all(math.isfinite(row["constraint_discriminant"]) for row in rows)
    assert min(row["constraint_discriminant"] for row in rows) > 0.0
    max_relative_kinetic_numerator = max(
        abs(row["kinetic_numerator"]) / max(row["kinetic_cancellation_scale"], 1.0e-300)
        for row in rows
    )
    assert max_relative_kinetic_numerator < 1.0e-11
    print("RFG-R extended-EFT scalar audit reproduced")
    print(f"grid points                         = {len(rows)}")
    print(f"max relative kinetic numerator      = {max_relative_kinetic_numerator:.3e}")
    print("result                              = degenerate scalar kinetic sector")


if __name__ == "__main__":
    main()
