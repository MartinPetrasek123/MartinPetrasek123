#!/usr/bin/env python3
"""Cross-check action identities without using the reconstruction residuals.

The main validator checks the constructed density and pressure.  This script
independently differentiates the declared background and reconstructs the
KGB alpha functions directly from the returned covariant action coefficients.
It deliberately avoids centered differences at the two smooth-future joins.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from ru_kgb import RUKGBParams, background


ROOT = Path(__file__).resolve().parents[1]


def relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(1.0, abs(left), abs(right))


def main() -> None:
    params = RUKGBParams()
    step = 2.0e-6
    # Do not straddle the C-infinity joins at a=1 and a=2.
    samples = (1.0e-6, 1.0e-4, 1.0e-2, 0.1, 0.5, 0.9, 1.1, 1.3, 1.5, 1.8, 1.9, 3.0, 10.0, 1.0e3)

    maximum_e_n = 0.0
    maximum_r_n = 0.0
    maximum_alpha_b = 0.0
    maximum_alpha_k = 0.0
    maximum_scalar_equation = 0.0

    for a in samples:
        row = background(a, params)
        n = math.log(a)
        plus = background(math.exp(n + step), params)
        minus = background(math.exp(n - step), params)

        e_n_finite_difference = (plus["E"] - minus["E"]) / (2.0 * step)
        r_n_finite_difference = (plus["rho_R"] - minus["rho_R"]) / (2.0 * step)

        e = row["E"]
        e_n = row["E_N"]
        kinetic = 0.5 * e * e
        a_hat = row["A_hat"]
        b_hat = row["B_hat"]
        c_hat = row["C_hat"]
        c_hat_phi = row["C_hat_phi"]
        potential = row["V_hat"]

        # For G2=A X+B X^2-V and G3=C X, these follow directly from
        # the covariant KGB energy-momentum tensor and alpha definitions.
        density = (
            a_hat * kinetic
            + 3.0 * b_hat * kinetic * kinetic
            + potential
            + 6.0 * e * e * kinetic * c_hat
            - 2.0 * c_hat_phi * kinetic * kinetic
        )
        pressure = (
            a_hat * kinetic
            + b_hat * kinetic * kinetic
            - potential
            - 2.0 * c_hat_phi * kinetic * kinetic
            - 2.0 * c_hat * kinetic * e * e_n
        )
        alpha_b = e * e * c_hat
        alpha_k = (
            2.0 * kinetic * a_hat
            + 12.0 * b_hat * kinetic * kinetic
            - 8.0 * c_hat_phi * kinetic * kinetic
            + 12.0 * e * e * kinetic * c_hat
        ) / (e * e)

        maximum_e_n = max(maximum_e_n, relative_difference(e_n_finite_difference, e_n))
        maximum_r_n = max(maximum_r_n, relative_difference(r_n_finite_difference, row["rho_R_N"]))
        maximum_alpha_b = max(maximum_alpha_b, relative_difference(alpha_b, row["alpha_B"]))
        maximum_alpha_k = max(maximum_alpha_k, relative_difference(alpha_k, row["alpha_K"]))

        # At a < 1e-2 the exact density is a cancellation of O(a^-4) terms.
        # The primary validator checks that regime at 80 digits.  Here the
        # double-precision finite-difference form is used only where stable.
        if a >= 1.0e-2:
            scalar_equation = (
                (plus["rho_phi_reconstructed"] - minus["rho_phi_reconstructed"]) / (2.0 * step)
                + 3.0 * (density + pressure)
            )
            maximum_scalar_equation = max(
                maximum_scalar_equation,
                abs(scalar_equation) / max(1.0, abs(density), abs(pressure)),
            )

    assert maximum_e_n < 2.0e-8
    assert maximum_r_n < 2.0e-8
    assert maximum_alpha_b < 2.0e-13
    assert maximum_alpha_k < 2.0e-12
    assert maximum_scalar_equation < 2.0e-8

    report = {
        "scope": (
            "Independent finite-difference and action-identity check for the declared "
            "G2+G3 KGB action. It is not a CMB, matter, or local-gravity likelihood."
        ),
        "samples": list(samples),
        "max_relative_E_N_finite_difference_residual": maximum_e_n,
        "max_relative_rho_R_N_finite_difference_residual": maximum_r_n,
        "max_relative_alpha_B_action_residual": maximum_alpha_b,
        "max_relative_alpha_K_action_residual": maximum_alpha_k,
        "max_relative_scalar_equation_residual": maximum_scalar_equation,
    }
    output = ROOT / "generated" / "independent_action_identities.json"
    output.write_text(json.dumps(report, indent=2) + "\n")

    print("Independent KGB action-identity validation OK")
    for key, value in report.items():
        if key.startswith("max_"):
            print(f"{key} = {value:.3e}")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
