#!/usr/bin/env python3
"""Quasi-static matter-growth prediction of the explicit R-alpha KGB action.

This is the subhorizon, linear-matter sector of the stated action.  It is not
used as a substitute for a photon-baryon Boltzmann calculation.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from ru_kgb import RUKGBParams, background


ROOT = Path(__file__).resolve().parents[1]


def rhs(n: float, state: np.ndarray, params: RUKGBParams, kgb: bool) -> np.ndarray:
    a = float(np.exp(n))
    if kgb:
        row = background(a, params)
        e_n_over_e = row["E_N"] / row["E"]
        omega_m = row["Omega_m"]
        # For alpha_M=alpha_T=0, the luminal Horndeski high-k response is
        # mu_infinity=1+alpha_B^2/(2 D c_s^2), with eta_infinity=1.
        mu = 1.0 + row["alpha_B"] ** 2 / (2.0 * row["D"] * row["c_s2"])
    else:
        matter = params.omega_m0 * a**-3
        radiation = params.omega_r0 * a**-4
        e2 = matter + radiation + params.omega_R0
        e_n_over_e = (-3.0 * matter - 4.0 * radiation) / (2.0 * e2)
        omega_m = matter / e2
        mu = 1.0
    d, d_n = state
    return np.array([d_n, -(2.0 + e_n_over_e) * d_n + 1.5 * omega_m * mu * d])


def integrate(params: RUKGBParams, kgb: bool, a_initial: float = 1.0e-3, steps: int = 3000) -> tuple[np.ndarray, np.ndarray]:
    n = np.linspace(np.log(a_initial), 0.0, steps + 1)
    h = n[1] - n[0]
    state = np.array([a_initial, a_initial])
    out = np.empty((steps + 1, 2))
    out[0] = state
    for i in range(steps):
        k1 = rhs(n[i], state, params, kgb)
        k2 = rhs(n[i] + 0.5 * h, state + 0.5 * h * k1, params, kgb)
        k3 = rhs(n[i] + 0.5 * h, state + 0.5 * h * k2, params, kgb)
        k4 = rhs(n[i] + h, state + h * k3, params, kgb)
        state = state + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        out[i + 1] = state
    return n, out


def main() -> None:
    params = RUKGBParams()
    n, kgb_solution = integrate(params, kgb=True)
    _, lcdm_solution = integrate(params, kgb=False)
    d_kgb = kgb_solution[:, 0] / kgb_solution[-1, 0]
    d_lcdm = lcdm_solution[:, 0] / lcdm_solution[-1, 0]
    f_kgb = kgb_solution[:, 1] / kgb_solution[:, 0]
    f_lcdm = lcdm_solution[:, 1] / lcdm_solution[:, 0]

    rows = []
    for ni, dk, dl, fk, fl in zip(n, d_kgb, d_lcdm, f_kgb, f_lcdm):
        a = float(np.exp(ni))
        row = background(a, params)
        mu = 1.0 + row["alpha_B"] ** 2 / (2.0 * row["D"] * row["c_s2"])
        rows.append({
            "a": a,
            "z": 1.0 / a - 1.0,
            "D_kgb_normalized": float(dk),
            "D_lcdm_normalized": float(dl),
            "D_ratio": float(dk / dl),
            "f_kgb": float(fk),
            "f_lcdm": float(fl),
            "fD_kgb": float(fk * dk),
            "fD_lcdm": float(fl * dl),
            "mu_infinity": float(mu),
            "eta_infinity": 1.0,
            "Sigma_infinity": float(mu),
        })

    output = ROOT / "generated" / "matter_qs_prediction.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    max_mu = max(row["mu_infinity"] for row in rows)
    min_ratio = min(row["D_ratio"] for row in rows)
    max_ratio = max(row["D_ratio"] for row in rows)
    summary = (
        "Quasi-static KGB growth prediction: "
        f"max mu_infinity={max_mu:.8f}, D/D_LCDM in [{min_ratio:.8f}, {max_ratio:.8f}]."
    )
    (ROOT / "generated" / "matter_qs_summary.txt").write_text(summary + "\n")
    print(summary)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
