#!/usr/bin/env python3
"""Exact ADM-to-extended-EFT map for the RFG-R cosmological action.

The map follows Eqs. (7)--(19) of Frusciante, Papadomanolakis and
Silvestri, arXiv:1601.04064. Their extrinsic curvature has the opposite sign
to the RFG-R convention, so the calculation first applies
``K_reference=-K_RFG``. It deliberately retains the extended operator
``bar_m5 delta R^(3) delta K``. Dropping that coefficient would change the
linear scalar theory and would not be an RFG-R calculation.

Units in the generated table are powers of the RFG-R H0 and M_Pl:

    c_hat          = c / (M_Pl^2 H0^2)
    Lambda_hat     = Lambda / (M_Pl^2 H0^2)
    M2_4_hat       = M2^4 / (M_Pl^2 H0^2)
    M1_3_hat       = bar_M1^3 / (M_Pl^2 H0)
    M2_bar_hat     = bar_M2^2 / M_Pl^2
    m5_bar_hat     = H0 bar_m5 / M_Pl^2

The table deliberately records basis-level EFT coefficients rather than
solver-specific gamma normalizations. The separate m5_bar_hat column is
essential. A backend which cannot consume it cannot produce an exact RFG-R
CMB or matter spectrum.
"""

from __future__ import annotations

import csv
from functools import lru_cache
import math
from pathlib import Path

import numpy as np

from rfg_regularized import (
    RFGRegularizedParams,
    Q,
    Q_prime,
    Q_second,
    Q_third,
    potential,
    potential_prime,
    potential_second,
    potential_third,
    response_prime,
    response_second,
    solve_E,
)


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "generated" / "tables"


def background_derivatives(a: float, params: RFGRegularizedParams) -> dict[str, float]:
    """Return E and exact first/two-log-a derivatives of the implicit branch."""
    x = solve_E(a, params)
    matter = params.omega_m0 * a**-3
    radiation = params.omega_r0 * a**-4
    source_n = -3.0 * matter - 4.0 * radiation
    source_nn = 9.0 * matter + 16.0 * radiation
    denominator = 2.0 * x - response_prime(x, params)
    x_n = source_n / denominator
    x_nn = (source_nn - (2.0 - response_second(x, params)) * x_n * x_n) / denominator
    h_dot = x * x_n
    h_ddot = x * (x_n * x_n + x * x_nn)
    return {
        "E": x,
        "dE_dln_a": x_n,
        "d2E_dln_a2": x_nn,
        "Hdot_over_H0_sq": h_dot,
        "Hddot_over_H0_cu": h_ddot,
        "background_denominator": denominator,
    }


@lru_cache(maxsize=None)
def extended_eft_coefficients(a: float, params: RFGRegularizedParams) -> dict[str, float]:
    """Map the RFG-R ADM action to all nonzero extended-EFT coefficients.

    The starting ADM Lagrangian is

        L = M_Pl^2/2 [Q(X)(R3 + S - K^2) + 2 H0^2 V(X)], X=-K/(3H0).

    It has L_N=L_NN=L_NR=L_NS=0. The C coefficient in the general ADM map
    is -M_Pl^2 Q_X/(6 H0), so bar_m5 is nonzero whenever Q_X is nonzero.
    """
    d = background_derivatives(a, params)
    x = d["E"]
    q = Q(x, params)
    q_x = Q_prime(x, params)
    q_xx = Q_second(x, params)
    v = potential(x, params)
    v_x = potential_prime(x, params)
    v_xx = potential_second(x, params)
    h_dot = d["Hdot_over_H0_sq"]
    h_ddot = d["Hddot_over_H0_cu"]

    # F = L_K - 2 H L_S, in units M_Pl^2 H0. K here is the
    # reference convention, K=-K_RFG, hence X=-K/(3 H0).
    f_hat = x * x * q_x + 2.0 * x * q - v_x / 3.0
    q_minus_one = q - 1.0
    # Keep the GR piece separate. At early times it cancels identically from
    # c and Lambda, while evaluating the cancellation in IEEE double precision
    # would erase the small RFG-R correction.
    f_mod_hat = x * x * q_x + 2.0 * x * q_minus_one - v_x / 3.0
    f_mod_x_hat = 2.0 * q_minus_one + 4.0 * x * q_x + x * x * q_xx - v_xx / 3.0
    f_mod_dot_hat = f_mod_x_hat * h_dot
    f_dot_hat = 2.0 * h_dot + f_mod_dot_hat
    q_dot_hat = q_x * h_dot
    q_ddot_hat = q_xx * h_dot * h_dot + q_x * h_ddot

    # Eq. (19) of arXiv:1601.04064, with H0=M_Pl=1 during evaluation.
    c_hat = 0.5 * f_mod_dot_hat + 0.5 * x * q_dot_hat - 0.5 * q_ddot_hat - q_minus_one * h_dot
    lambda_hat = (
        v
        + f_mod_dot_hat
        + 3.0 * x * f_mod_hat
        - 6.0 * q_minus_one * x * x
        - q_ddot_hat
        - 2.0 * x * q_dot_hat
        - 2.0 * q_minus_one * h_dot
    )
    m2_4_hat = -0.5 * c_hat
    m1_3_hat = -q_dot_hat
    m2_bar_hat = (x * x * q_xx) / 3.0 + (4.0 * x * q_x) / 3.0 - v_xx / 9.0
    m3_bar_hat = 0.0
    m5_bar_hat = -q_x / 3.0

    return {
        "a": a,
        "z": 1.0 / a - 1.0,
        **d,
        "Q": q,
        "Q_X": q_x,
        "Q_XX": q_xx,
        "V": v,
        "V_X": v_x,
        "V_XX": v_xx,
        "Omega": q - 1.0,
        "c_hat": c_hat,
        "Lambda_hat": lambda_hat,
        "M2_4_hat": m2_4_hat,
        "M1_3_hat": m1_3_hat,
        "M2_bar_hat": m2_bar_hat,
        "M3_bar_hat": m3_bar_hat,
        "Mhat2_hat": 0.0,
        "m2_sq_hat": 0.0,
        "m5_bar_hat": m5_bar_hat,
        "deltaR_deltaK_coefficient_hat": 0.5 * m5_bar_hat,
        "F_hat": f_hat,
        "Fdot_hat": f_dot_hat,
    }


@lru_cache(maxsize=None)
def extended_eft_w_coefficients(a: float, params: RFGRegularizedParams) -> dict[str, float]:
    """Return W_i and their analytic cosmic-time derivatives.

    The simplifications below use the RFG-R identities
    ``M2_4=-c/2`` and ``M1_3=-dot(Q)`` before evaluation, avoiding any
    cancellation between independently rounded EFT entries.  No finite
    difference in time or scale factor is used.
    """
    row = extended_eft_coefficients(a, params)
    h = float(row["E"])
    h_dot = float(row["Hdot_over_H0_sq"])
    h_ddot = float(row["Hddot_over_H0_cu"])
    q = float(row["Q"])
    q_x = float(row["Q_X"])
    q_xx = float(row["Q_XX"])
    q_xxx = Q_third(h, params)
    v_xx = float(row["V_XX"])
    v_xxx = potential_third(h, params)
    q_dot = q_x * h_dot
    m2 = float(row["M2_bar_hat"])
    m2_dot = (2.0 * h * q_xx + h * h * q_xxx + 4.0 * q_x + 4.0 * h * q_xx) * h_dot / 3.0 - v_xxx * h_dot / 9.0
    m5 = -q_x / 3.0
    m5_dot = -q_xx * h_dot / 3.0
    m5_ddot = -(q_xxx * h_dot * h_dot + q_xx * h_ddot) / 3.0

    w0_numerator = q + 3.0 * h * m5 + 3.0 * m5_dot
    w0_numerator_dot = q_dot + 3.0 * (h_dot * m5 + h * m5_dot) + 3.0 * m5_ddot
    w0 = -w0_numerator / (a * a)
    w0_dot = -(w0_numerator_dot - 2.0 * h * w0_numerator) / (a * a)

    w1 = -3.0 * h * h * q - 4.5 * h * h * m2
    w1_dot = -6.0 * h * h_dot * q - 3.0 * h * h * q_dot - 9.0 * h * h_dot * m2 - 4.5 * h * h * m2_dot

    w4_numerator = -2.0 * h * q - 3.0 * h * m2
    w4_numerator_dot = -2.0 * h_dot * q - 2.0 * h * q_dot - 3.0 * h_dot * m2 - 3.0 * h * m2_dot
    w4 = w4_numerator / (a * a)
    w4_dot = (w4_numerator_dot - 2.0 * h * w4_numerator) / (a * a)

    w5_numerator = 2.0 * q + 3.0 * m2
    w5_numerator_dot = 2.0 * q_dot + 3.0 * m2_dot
    w5 = w5_numerator / (a * a)
    w5_dot = (w5_numerator_dot - 2.0 * h * w5_numerator) / (a * a)

    w6_numerator = -2.0 * q - 6.0 * h * m5
    w6_numerator_dot = -2.0 * q_dot - 6.0 * (h_dot * m5 + h * m5_dot)
    w6 = w6_numerator / (a * a)
    w6_dot = (w6_numerator_dot - 2.0 * h * w6_numerator) / (a * a)

    w7 = -m2 / (2.0 * a**4)
    w7_dot = -(m2_dot - 4.0 * h * m2) / (2.0 * a**4)
    return {
        **row,
        "Q_XXX": q_xxx,
        "V_XXX": v_xxx,
        "Qdot": q_dot,
        "M2_bar_dot_hat": m2_dot,
        "m5_bar_dot_hat": m5_dot,
        "W0": w0,
        "W1": w1,
        "W4": w4,
        "W5": w5,
        "W6": w6,
        "W7": w7,
        "W0_dot": w0_dot,
        "W1_dot": w1_dot,
        "W4_dot": w4_dot,
        "W5_dot": w5_dot,
        "W6_dot": w6_dot,
        "W7_dot": w7_dot,
    }


def _write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_summary(path: Path, rows: list[dict[str, float]], params: RFGRegularizedParams) -> None:
    m5_values = [row["m5_bar_hat"] for row in rows]
    q_values = [row["Q"] for row in rows]
    text = [
        "# RFG-R Extended EFT Map",
        "",
        "The table is the exact ADM-to-EFT map of the RFG-R cosmological action.",
        "It retains the extended operator bar_m5 deltaR3 deltaK.",
        "",
        f"theta = {params.theta:.8f}",
        f"epsilon = {params.epsilon:.1e}",
        f"p = {params.p}",
        f"a range = [{rows[0]['a']:.1e}, {rows[-1]['a']:.1e}]",
        f"min Q = {min(q_values):.12e}",
        f"max Q = {max(q_values):.12e}",
        f"min m5_bar_hat = {min(m5_values):.12e}",
        f"max m5_bar_hat = {max(m5_values):.12e}",
        "",
        "m5_bar_hat = H0 bar_m5 / M_Pl^2 = -Q_X/3.",
        "This coefficient is nonzero on the reference cosmological branch.",
        "Stock H-EFTCAMB exposes gamma_1...gamma_6 but no bar_m5 deltaR3 deltaK",
        "operator, so the table is an exact theory input rather than a completed",
        "stock-H-EFTCAMB spectrum run.",
        "The companion pure-gravity scalar audit is degenerate on its 2,401-point",
        "(a,k) grid; a sourced multi-fluid reduction is required before spectra",
        "or CMB/matter likelihoods can be evaluated.",
    ]
    path.write_text("\n".join(text) + "\n", encoding="utf-8")


def main() -> None:
    params = RFGRegularizedParams()
    params.validate()
    TABLES.mkdir(parents=True, exist_ok=True)
    rows = [extended_eft_coefficients(float(a), params) for a in np.logspace(-7, 0, 81)]
    if not all(math.isfinite(value) for row in rows for value in row.values()):
        raise RuntimeError("non-finite value in extended EFT map")
    _write_csv(TABLES / "extended_eft_mapping.csv", rows)
    _write_summary(TABLES / "extended_eft_mapping_summary.md", rows, params)
    print(f"Generated {TABLES / 'extended_eft_mapping.csv'}")
    print(f"Generated {TABLES / 'extended_eft_mapping_summary.md'}")


if __name__ == "__main__":
    main()
