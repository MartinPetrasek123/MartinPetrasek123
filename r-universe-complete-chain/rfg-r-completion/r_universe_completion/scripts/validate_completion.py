#!/usr/bin/env python3
"""Independent mathematical and numerical checks for the RFG-R completion."""

from __future__ import annotations

import math

import numpy as np

from rfg_regularized import (
    RFGRegularizedParams,
    Q,
    background_residual,
    observables,
    original_Q,
    original_response,
    potential,
    potential_small_x_coefficient,
    response,
    source_for_potential,
)
from ppn_likelihood import W_LOCAL_GR, cosmology_weight, ppn_parameters, solar_weyl_ratio


def main() -> None:
    p = RFGRegularizedParams()
    p.validate()

    # Exact local Einstein-Hilbert coefficients.
    assert Q(0.0, p) == 1.0
    assert potential(0.0, p) == 0.0

    # Recovery of the original branch above the cosmological lower bound X=0.8.
    xs = np.logspace(math.log10(0.8), 5.0, 300)
    response_errors = [
        abs(response(float(x), p) / original_response(float(x), p) - 1.0) for x in xs
    ]
    q_errors = [
        abs((Q(float(x), p) - original_Q(float(x), p)) / (original_Q(float(x), p) - 1.0))
        for x in xs
    ]
    max_response_error = max(response_errors)
    max_q_relative_error = max(q_errors)
    # The analytic correction is O((epsilon/X)^p) <= 2.5e-32 here; double
    # precision arithmetic itself saturates near 1e-15 at the largest X.
    assert max_response_error < 1.0e-12
    assert max_q_relative_error < 1.0e-12

    # The reconstructed V satisfies V - X V_X = 3F.  A centered numerical derivative
    # is kept deliberately independent of the identity used in potential_prime().
    residuals = []
    for x in [1.0e-5, 1.0e-3, 0.1, 0.8, 1.0, 10.0]:
        h = max(1.0e-5 * x, 1.0e-9)
        vp = (potential(x + h, p) - potential(x - h, p)) / (2.0 * h)
        residuals.append(abs(potential(x, p) - x * vp - 3.0 * source_for_potential(x, p)))
    max_potential_residual = max(residuals)
    assert max_potential_residual < 2.0e-5

    # The predicted small-X leading power is X^(p+2), so it cannot alter the GR
    # quadratic ADM kinetic operator. Compare at a point safely above roundoff.
    x_small = p.epsilon * 1.0e-3
    low_x_ratio = potential(x_small, p) / (potential_small_x_coefficient(p) * x_small ** (p.p + 2))
    assert abs(low_x_ratio - 1.0) < 2.0e-3

    # Cosmological branch closure and positive tensor normalization.
    rows = [observables(float(a), p) for a in np.logspace(-8, 2, 280)]
    max_branch_residual = max(
        abs(background_residual(r["E"], r["a"], p)) / (r["E"] * r["E"]) for r in rows
    )
    max_closure_error = max(abs(r["closure"] - 1.0) for r in rows)
    min_qt = min(r["Q_T"] for r in rows)
    assert max_branch_residual < 1.0e-12
    assert max_closure_error < 1.0e-11
    assert min_qt > 0.0

    solar_ratio = solar_weyl_ratio(149597870700.0)
    assert solar_ratio > W_LOCAL_GR
    assert cosmology_weight(solar_ratio) == 0.0
    local_ppn = ppn_parameters(solar_ratio)
    assert local_ppn == {"gamma": 1.0, "beta": 1.0, "alpha1": 0.0, "alpha2": 0.0}

    print("RFG-R completion validation OK")
    print(f"max high-X response error  = {max_response_error:.3e}")
    print(f"max high-X Q error         = {max_q_relative_error:.3e}")
    print(f"max potential ODE residual = {max_potential_residual:.3e}")
    print(f"low-X potential ratio      = {low_x_ratio:.10f}")
    print(f"max branch residual        = {max_branch_residual:.3e}")
    print(f"max density closure error  = {max_closure_error:.3e}")
    print(f"min Q_T                    = {min_qt:.10f}")
    print(f"solar Weyl ratio at 1 AU    = {solar_ratio:.3e}")
    print("PPN result in GR domain    = gamma=1, beta=1, alpha1=0, alpha2=0")


if __name__ == "__main__":
    main()
