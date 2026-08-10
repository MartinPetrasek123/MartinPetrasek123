#!/usr/bin/env python3
"""Scalar-sector audit for the exact RFG-R extended-EFT action.

This evaluates the reduced *pure-gravity* unitary-gauge action of
Frusciante, Papadomanolakis and Silvestri, arXiv:1601.04064,
Eqs. (85)--(86), for the complete RFG-R ADM-to-EFT map.  In particular it
keeps the nonzero ``bar_m5 delta R^(3) delta K`` term; setting it to zero is
not an RFG-R calculation.

The result is a gravity-sector diagnostic, not a CMB or matter likelihood.
The cited quadratic action is derived after omitting matter perturbations. A
degenerate result blocks a standalone EFT/Boltzmann evolution; a full
multi-species reduction must determine whether sourced matter constraints lift
the degeneracy or impose an additional physical constraint.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from extended_eft_mapping import extended_eft_coefficients
from rfg_regularized import RFGRegularizedParams


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "generated" / "tables"


def _w_coefficients(a: float, params: RFGRegularizedParams) -> dict[str, float]:
    """Return W_i of Eqs. (85)--(86), in units H0=M_Pl=1."""
    row = extended_eft_coefficients(a, params)
    h = row["E"]
    q = row["Q"]
    omega_dot = row["Q_X"] * row["Hdot_over_H0_sq"]
    # bar_m5=-Q_X/3, therefore dot(bar_m5)=-Q_XX dot(H)/3 exactly.
    # This is a defining action derivative, not a finite-difference input.
    m5_dot = -row["Q_XX"] * row["Hdot_over_H0_sq"] / 3.0

    w0 = -(q + 3.0 * h * row["m5_bar_hat"] + 3.0 * m5_dot) / a**2
    w1 = (
        row["c_hat"]
        + 2.0 * row["M2_4_hat"]
        - 3.0 * h * h * q
        - 3.0 * h * omega_dot
        - 1.5 * h * h * row["M3_bar_hat"]
        - 4.5 * h * h * row["M2_bar_hat"]
        - 3.0 * h * row["M1_3_hat"]
    )
    w4 = (
        -2.0 * h * q
        - omega_dot
        - h * row["M3_bar_hat"]
        - row["M1_3_hat"]
        - 3.0 * h * row["M2_bar_hat"]
    ) / a**2
    w5 = (2.0 * q + row["M3_bar_hat"] + 3.0 * row["M2_bar_hat"]) / a**2
    w6 = (-2.0 * q - 6.0 * h * row["m5_bar_hat"]) / a**2
    w7 = -(row["M3_bar_hat"] + row["M2_bar_hat"]) / (2.0 * a**4)
    return {"W0": w0, "W1": w1, "W4": w4, "W5": w5, "W6": w6, "W7": w7, **row}


def _reduced_coefficients(a: float, k_over_h0: float, params: RFGRegularizedParams) -> dict[str, float]:
    """Reduce Eq. (85) with the RFG-R nonzero EFT functions.

    RFG-R has m2^2=lambda_i=0 but has both bar_m5 and W7 nonzero.  Equations
    (110)--(114) therefore reduce to the expressions below.  ``k_over_h0``
    is the physical comoving wavenumber in H0 units.
    """
    values = _w_coefficients(a, params)
    w0, w1 = values["W0"], values["W1"]
    w4, w5, w6, w7 = values["W4"], values["W5"], values["W6"], values["W7"]
    m5 = values["m5_bar_hat"]
    discriminant = w4 * w4 - 4.0 * w1 * w7
    denominator = a * a * discriminant
    if denominator == 0.0:
        raise ZeroDivisionError("singular lapse/shift constraint")

    kinetic_prefactor = 6.0 * a * a * w7 + w5
    kinetic_term_1 = 3.0 * a**4 * w4 * w4
    kinetic_term_2 = 2.0 * a * a * w1 * w5
    kinetic_numerator = kinetic_prefactor * (kinetic_term_1 + kinetic_term_2)
    kinetic = kinetic_numerator / (2.0 * denominator)
    b_bar = (
        a * a * w0 * discriminant
        + k_over_h0**2 / a**6 * (-2.0 * a**4 * m5 * w4 * w6 - 4.0 * m5 * m5 * w1)
    ) / denominator
    v_bar = -(
        6.0 * a**4 * w4 * w7 * w6
        + a * a * w4 * w5 * w6
        + 6.0 * m5 * w4 * w4
        + 4.0 * m5 * w1 * w5 / a**2
    ) / denominator
    cancellation_scale = abs(kinetic_prefactor) * (abs(kinetic_term_1) + abs(kinetic_term_2))
    scalar_kinetic_degenerate = abs(kinetic_numerator) <= 1.0e-11 * max(cancellation_scale, 1.0e-300)
    return {
        "kinetic": kinetic,
        "B_bar": b_bar,
        "V_bar": v_bar,
        "constraint_discriminant": discriminant,
        "kinetic_numerator": kinetic_numerator,
        "kinetic_cancellation_scale": cancellation_scale,
        "scalar_kinetic_degenerate": scalar_kinetic_degenerate,
        **values,
    }


def _gradient(a: float, k_over_h0: float, params: RFGRegularizedParams) -> float:
    """Reject an undefined pure-gravity sound-speed calculation.

    The audited RFG-R branch has an identically degenerate scalar kinetic
    coefficient. Evaluating a putative gradient by numerically differentiating
    another reduced coefficient would not create a physical propagation speed.
    A nondegenerate action would require its own analytic derivative map.
    """
    del a, k_over_h0, params
    raise RuntimeError("the pure-gravity scalar gradient is undefined on the degenerate RFG-R branch")


def scalar_stability_row(a: float, k_over_h0: float, params: RFGRegularizedParams) -> dict[str, float]:
    """Return the no-ghost and scalar-gradient diagnostics at one (a,k)."""
    reduced = _reduced_coefficients(a, k_over_h0, params)
    kinetic = reduced["kinetic"]
    if reduced["scalar_kinetic_degenerate"]:
        # The reduced quadratic action has no zeta-dot squared term.  A sound
        # speed is therefore undefined, rather than very large or very small.
        gradient = float("nan")
        sound_speed_sq = float("nan")
    else:
        gradient = _gradient(a, k_over_h0, params)
        sound_speed_sq = gradient / kinetic
    return {
        "a": a,
        "z": 1.0 / a - 1.0,
        "k_over_H0": k_over_h0,
        "kinetic_L_zeta_dot_zeta": kinetic,
        "gradient_G": gradient,
        "sound_speed_sq": sound_speed_sq,
        "constraint_discriminant": reduced["constraint_discriminant"],
        "kinetic_numerator": reduced["kinetic_numerator"],
        "kinetic_cancellation_scale": reduced["kinetic_cancellation_scale"],
        "scalar_kinetic_degenerate": int(reduced["scalar_kinetic_degenerate"]),
        "Q": reduced["Q"],
        "m5_bar_hat": reduced["m5_bar_hat"],
    }


def main() -> None:
    params = RFGRegularizedParams()
    params.validate()
    TABLES.mkdir(parents=True, exist_ok=True)
    scale_factors = np.logspace(-7.0, 0.0, 49)
    wavenumbers = np.logspace(-4.0, 5.0, 49)
    rows = [scalar_stability_row(float(a), float(k), params) for a in scale_factors for k in wavenumbers]
    finite_fields = ("a", "z", "k_over_H0", "kinetic_L_zeta_dot_zeta", "constraint_discriminant", "Q", "m5_bar_hat")
    if not all(math.isfinite(row[field]) for row in rows for field in finite_fields):
        raise RuntimeError("non-finite background or constraint value in scalar stability audit")

    path = TABLES / "extended_eft_scalar_stability.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    min_kinetic = min(row["kinetic_L_zeta_dot_zeta"] for row in rows)
    min_constraint = min(row["constraint_discriminant"] for row in rows)
    degenerate_rows = [row for row in rows if row["scalar_kinetic_degenerate"]]
    print(f"Generated {path}")
    print(f"min L_zeta_dot_zeta = {min_kinetic:.12e}")
    print(f"min constraint disc  = {min_constraint:.12e}")
    print(f"kinetic-degenerate points = {len(degenerate_rows)} / {len(rows)}")

    if min_constraint <= 0.0:
        raise RuntimeError("singular/negative lapse-shift constraint discriminant")
    if degenerate_rows:
        raise RuntimeError(
            "RFG-R has a degenerate pure-gravity quadratic scalar kinetic term on the audited branch; "
            "the sourced multi-fluid reduction and kinetic hierarchy must be implemented in a full solver before a CMB/matter likelihood."
        )
    if min_kinetic <= 0.0:
        raise RuntimeError("scalar ghost found in gravity-sector audit")
    min_gradient = min(row["gradient_G"] for row in rows)
    min_sound_speed_sq = min(row["sound_speed_sq"] for row in rows)
    print(f"min G               = {min_gradient:.12e}")
    print(f"min c_s^2           = {min_sound_speed_sq:.12e}")
    if min_sound_speed_sq <= 0.0:
        raise RuntimeError("scalar gradient instability found in gravity-sector audit")


if __name__ == "__main__":
    main()
