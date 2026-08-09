#!/usr/bin/env python3
"""Generate RFG-R tables and figures used by the standalone paper."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from rfg_regularized import (
    RFGRegularizedParams,
    Q,
    eft_coefficients,
    find_z_acc,
    observables,
    original_Q,
    original_response,
    response,
)


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "generated" / "tables"
FIGURES = ROOT / "generated" / "figures"


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(p: RFGRegularizedParams, rows: list[dict[str, float]]) -> None:
    now = observables(1.0, p)
    z1 = observables(0.5, p)
    zacc = find_z_acc(p)
    text = [
        "# RFG-R Reference Point",
        "",
        f"Omega_m0 = {p.omega_m0:.8f}",
        f"Omega_r0 = {p.omega_r0:.8f}",
        f"Omega_R0 = {p.omega_R0:.8f}",
        f"theta = {p.theta:.8f}",
        f"epsilon = {p.epsilon:.1e}",
        f"p = {p.p}",
        "",
        f"Q_T(0) = {now['Q_T']:.10f}",
        f"alpha_M(0) = {now['alpha_M']:.10f}",
        f"dL_GW/dL_EM(z=1) = {z1['dL_GW_over_dL_EM']:.10f}",
        f"z_acc = {zacc:.10f}",
        "",
        "The RFG-R values equal the original RFG values to numerical precision on the",
        "cosmological branch. Epsilon only resolves the otherwise singular X=0 limit.",
    ]
    (TABLES / "reference_point.md").write_text("\n".join(text) + "\n", encoding="utf-8")


def plot_regularization(p: RFGRegularizedParams) -> None:
    x = np.logspace(-12, 5, 700)
    ratio_response = np.array([response(float(v), p) / original_response(float(v), p) for v in x])
    q_change = np.array([abs(Q(float(v), p) - original_Q(float(v), p)) for v in x])
    plt.figure(figsize=(7.2, 4.6))
    plt.semilogx(x, ratio_response, label="R_epsilon / R_original")
    plt.semilogx(x, q_change, label="|Q_epsilon - Q_original|")
    plt.axvline(p.epsilon, color="black", linestyle="--", linewidth=0.9, label="epsilon")
    plt.ylim(-0.02, 1.05)
    plt.xlabel("X")
    plt.ylabel("dimensionless response")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGURES / "regularization_recovery.png", dpi=180)
    plt.close()


def plot_eft(rows: list[dict[str, float]]) -> None:
    a = np.array([r["a"] for r in rows])
    qt = np.array([r["Q_T"] for r in rows])
    am = np.array([r["alpha_M"] for r in rows])
    plt.figure(figsize=(7.2, 4.6))
    plt.semilogx(a, qt, label="Q_T")
    plt.semilogx(a, am, label="alpha_M")
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.xlabel("a")
    plt.ylabel("value")
    plt.xlim(1.0e-8, 1.0)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGURES / "tensor_and_planck_running.png", dpi=180)
    plt.close()


def plot_background(rows: list[dict[str, float]]) -> None:
    a = np.array([r["a"] for r in rows])
    om = np.array([r["Omega_m"] for r in rows])
    orad = np.array([r["Omega_r"] for r in rows])
    ore = np.array([r["Omega_R"] for r in rows])
    plt.figure(figsize=(7.2, 4.6))
    plt.semilogx(a, om, label="Omega_m")
    plt.semilogx(a, orad, label="Omega_r")
    plt.semilogx(a, ore, label="Omega_R")
    plt.xlabel("a")
    plt.ylabel("fractional density")
    plt.xlim(1.0e-8, 1.0)
    plt.ylim(0.0, 1.02)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGURES / "fractional_densities.png", dpi=180)
    plt.close()


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    p = RFGRegularizedParams()
    a_values = np.logspace(-8, 0, 301)
    rows = [observables(float(a), p) for a in a_values]
    write_csv(TABLES / "completion_observables.csv", rows)
    eft_rows = [eft_coefficients(float(a), p) for a in np.logspace(-7, 0, 81)]
    write_csv(TABLES / "eft_coefficients.csv", eft_rows)
    write_summary(p, rows)
    plot_regularization(p)
    plot_eft(rows)
    plot_background(rows)
    print(f"Generated tables under {TABLES}")
    print(f"Generated figures under {FIGURES}")


if __name__ == "__main__":
    main()
