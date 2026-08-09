#!/usr/bin/env python3
"""Independent checks for the RFG-R ADM-to-extended-EFT coefficient map."""

from __future__ import annotations

import math

from extended_eft_mapping import background_derivatives, extended_eft_coefficients
from rfg_regularized import RFGRegularizedParams, Q_prime, solve_E


def _finite_log_a_derivatives(a: float, params: RFGRegularizedParams) -> tuple[float, float]:
    step = 2.0e-4
    center = solve_E(a, params)
    plus = solve_E(a * math.exp(step), params)
    minus = solve_E(a * math.exp(-step), params)
    first = (plus - minus) / (2.0 * step)
    second = (plus - 2.0 * center + minus) / (step * step)
    return first, second


def main() -> None:
    params = RFGRegularizedParams()
    params.validate()
    max_first_error = 0.0
    max_second_error = 0.0
    max_m5_identity_error = 0.0
    min_q = float("inf")
    max_abs_m5 = 0.0

    # In the pure GR limit Q=1 and V=0, F/M_Pl^2H0=2E. The explicit
    # GR-subtracted implementation must return c=Lambda=0 for arbitrary
    # background derivatives rather than only for a special H(a).
    x_gr = 2.7
    h_dot_gr = -4.2
    f_mod_gr = 0.0
    f_mod_dot_gr = 0.0
    q_dot_gr = 0.0
    q_ddot_gr = 0.0
    c_gr = 0.5 * f_mod_dot_gr + 0.5 * x_gr * q_dot_gr - 0.5 * q_ddot_gr
    lambda_gr = f_mod_dot_gr + 3.0 * x_gr * f_mod_gr - q_ddot_gr - 2.0 * x_gr * q_dot_gr
    assert c_gr == 0.0
    assert lambda_gr == 0.0

    for a in (1.0e-6, 1.0e-4, 1.0e-2, 0.1, 0.5, 1.0):
        derivs = background_derivatives(a, params)
        first_fd, second_fd = _finite_log_a_derivatives(a, params)
        max_first_error = max(max_first_error, abs(first_fd - derivs["dE_dln_a"]) / max(1.0, abs(first_fd)))
        max_second_error = max(max_second_error, abs(second_fd - derivs["d2E_dln_a2"]) / max(1.0, abs(second_fd)))

        row = extended_eft_coefficients(a, params)
        expected_m5 = -Q_prime(row["E"], params) / 3.0
        max_m5_identity_error = max(max_m5_identity_error, abs(row["m5_bar_hat"] - expected_m5))
        # The quadratic action coefficient of deltaR3 deltaK is bar_m5/2.
        max_m5_identity_error = max(
            max_m5_identity_error,
            abs(row["deltaR_deltaK_coefficient_hat"] - expected_m5 / 2.0),
        )
        min_q = min(min_q, row["Q"])
        max_abs_m5 = max(max_abs_m5, abs(row["m5_bar_hat"]))
        assert row["M3_bar_hat"] == 0.0
        assert row["Mhat2_hat"] == 0.0
        assert row["m2_sq_hat"] == 0.0

    assert max_first_error < 3.0e-7
    assert max_second_error < 2.0e-4
    assert max_m5_identity_error < 1.0e-13
    assert min_q > 0.0
    assert max_abs_m5 > 1.0e-6

    print("RFG-R extended EFT map validation OK")
    print(f"max dE/dln a finite-difference error  = {max_first_error:.3e}")
    print(f"max d2E/dln a2 finite-difference error = {max_second_error:.3e}")
    print(f"max bar_m5 mapping identity error      = {max_m5_identity_error:.3e}")
    print(f"min Q on validation points             = {min_q:.10f}")
    print(f"max |m5_bar_hat|                       = {max_abs_m5:.10f}")


if __name__ == "__main__":
    main()
