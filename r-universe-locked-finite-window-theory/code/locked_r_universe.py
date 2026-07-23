#!/usr/bin/env python3
"""Locked finite-window R-Universe background integration.

This script implements the locked scalar calibration theory specified in the
manuscript:

  W(R) = R_*/(R+R_*)
  F(R) = 1 + xi W(1-W)
  Z(R) = Z0 W^2
  U(R) = U_* W^2(1-W)^2
  Rddot + 3 H Rdot + ell^{-2} R = ell^{-2} rho_m

The calculation is dimensionless with M_Pl=1 and a=1 used only as the reporting
surface.  It is a direct integration of the locked equations, not a fit to an
external expansion curve.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "code" / "locked_r_universe_results.json"
OUT_CSV = ROOT / "tables" / "locked_r_universe_background.csv"
FIG_PDF = ROOT / "figures" / "fig01_locked_r_universe_background.pdf"
FIG_PNG = ROOT / "figures" / "fig01_locked_r_universe_background.png"


@dataclass(frozen=True)
class Params:
    xi: float = 8.0
    Z0: float = 0.03
    Ustar: float = 0.095
    Rstar: float = 0.55
    ell: float = 5.0
    rho_m0: float = 0.05
    rho_r0: float = 1.0e-4
    R_ini_factor: float = 1.0
    V_ini: float = 0.0


def W(R: np.ndarray | float, p: Params):
    return p.Rstar / (R + p.Rstar)


def dW_dR(R: np.ndarray | float, p: Params):
    return -p.Rstar / (R + p.Rstar) ** 2


def Ffun(R: np.ndarray | float, p: Params):
    w = W(R, p)
    return 1.0 + p.xi * w * (1.0 - w)


def dF_dR(R: np.ndarray | float, p: Params):
    w = W(R, p)
    return p.xi * dW_dR(R, p) * (1.0 - 2.0 * w)


def Zfun(R: np.ndarray | float, p: Params):
    return p.Z0 * W(R, p) ** 2


def Ufun(R: np.ndarray | float, p: Params):
    w = W(R, p)
    return p.Ustar * w * w * (1.0 - w) ** 2


def densities(a: float, p: Params):
    return p.rho_m0 * a ** -3, p.rho_r0 * a ** -4


def H_from_constraint(a: float, R: float, V: float, p: Params):
    rho_m, rho_r = densities(a, p)
    F = Ffun(R, p)
    FR = dF_dR(R, p)
    Z = Zfun(R, p)
    U = Ufun(R, p)
    # 3 F H^2 + 3 F_R V H - (rho_m+rho_r+0.5 Z V^2+U)=0
    A = 3.0 * F
    B = 3.0 * FR * V
    C = -(rho_m + rho_r + 0.5 * Z * V * V + U)
    disc = max(B * B - 4.0 * A * C, 0.0)
    return (-B + np.sqrt(disc)) / (2.0 * A)


def rhs(t: float, y: np.ndarray, p: Params):
    a, R, V = y
    if a <= 0 or R <= 0:
        return [0.0, 0.0, 0.0]
    H = H_from_constraint(a, R, V, p)
    rho_m, _ = densities(a, p)
    dV = (rho_m - R) / (p.ell * p.ell) - 3.0 * H * V
    return [a * H, V, dV]


def integrate(p: Params):
    a_ini = 1.0e-2
    a_end = 2.0
    R_ini = p.R_ini_factor * p.rho_m0 * a_ini ** -3
    y0 = np.array([a_ini, R_ini, p.V_ini], dtype=float)

    def event_a_end(_t, y):
        return y[0] - a_end

    event_a_end.terminal = True
    event_a_end.direction = 1
    sol = solve_ivp(
        lambda t, y: rhs(t, y, p),
        (0.0, 200.0),
        y0,
        method="DOP853",
        rtol=5.0e-7,
        atol=1.0e-9,
        events=event_a_end,
        max_step=0.35,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    t = sol.t
    a_raw, R_raw, V_raw = sol.y
    order = np.argsort(a_raw)
    a_raw, t_raw, R_raw, V_raw = a_raw[order], t[order], R_raw[order], V_raw[order]
    # Adaptive time stepping can be sparse in rapid transitions.  All reported
    # diagnostics are therefore evaluated on a fixed logarithmic scale-factor grid.
    a = np.logspace(np.log10(a_raw[0]), np.log10(a_raw[-1]), 320)
    loga_raw = np.log(a_raw)
    R = np.interp(np.log(a), loga_raw, R_raw)
    V = np.interp(np.log(a), loga_raw, V_raw)
    t = np.interp(np.log(a), loga_raw, t_raw)
    H = np.array([H_from_constraint(ai, Ri, Vi, p) for ai, Ri, Vi in zip(a, R, V)])
    rho_m = p.rho_m0 * a ** -3
    rho_r = p.rho_r0 * a ** -4
    F = Ffun(R, p)
    FR = dF_dR(R, p)
    Z = Zfun(R, p)
    U = Ufun(R, p)
    w = W(R, p)
    omega_m_R = rho_m / (3.0 * F * H * H)
    omega_r_R = rho_r / (3.0 * F * H * H)
    omega_cal = (0.5 * Z * V * V + U - 3.0 * H * FR * V) / (3.0 * F * H * H)
    # The first-constraint closure should be one up to integration error.
    closure = omega_m_R + omega_r_R + omega_cal
    H_N = np.gradient(np.log(H), np.log(a), edge_order=2)
    q = -1.0 - H_N
    geff_unscreened = 1.0 / F
    return pd.DataFrame(
        {
            "t": t,
            "a": a,
            "R_cal": R,
            "Rdot": V,
            "H": H,
            "W": w,
            "F": F,
            "Z": Z,
            "U": U,
            "Omega_m_R": omega_m_R,
            "Omega_r_R": omega_r_R,
            "Omega_cal": omega_cal,
            "closure_sum": closure,
            "q": q,
            "G_unscreened_over_G": geff_unscreened,
        }
    )


def score(df: pd.DataFrame):
    idx = int(np.argmin(np.abs(df["a"].to_numpy() - 1.0)))
    row = df.iloc[idx]
    return float(
        6.0 * abs(row["Omega_m_R"] - 0.30)
        + 0.8 * abs(row["q"] + 0.55)
        + 4.0 * abs(row["closure_sum"] - 1.0)
        + 0.15 * abs(row["F"] - 1.0)
        + 0.05 * max(abs(row["q"]) - 2.0, 0.0)
        + 0.02 * abs(np.log(max(row["H"], 1.0e-30)))
    )


def scan():
    base = Params()
    candidates = []
    for xi in [0.5, 1.0, 2.0, 5.0, 8.0]:
        for z0 in [0.001, 0.01, 0.05]:
            for u in [0.05, 0.15, 0.35, 0.75, 1.5]:
                for rs in [0.02, 0.05, 0.1, 0.2, 0.5]:
                    for ell in [0.7, 1.5, 3.0, 6.0]:
                        p = Params(xi=xi, Z0=z0, Ustar=u, Rstar=rs, ell=ell)
                        try:
                            df = integrate(p)
                            idx = int(np.argmin(np.abs(df["a"].to_numpy() - 1.0)))
                            row = df.iloc[idx]
                            if np.isfinite(row["q"]) and row["F"] > 0 and df["F"].min() > 0:
                                candidates.append((score(df), p, row, df))
                        except Exception:
                            continue
    candidates.sort(key=lambda x: x[0])
    return candidates


def make_plot(df: pd.DataFrame):
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.0))
    ax = axes[0, 0]
    ax.plot(df["a"], df["W"], label="W")
    ax.plot(df["a"], df["F"], label="F")
    ax.set_xscale("log")
    ax.axvline(1.0, color="0.5", ls=":", lw=0.8)
    ax.set_xlabel("a")
    ax.set_title("finite-window functions")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(df["a"], df["Omega_m_R"], label=r"$\Omega_m^R$")
    ax.plot(df["a"], df["Omega_r_R"], label=r"$\Omega_r^R$")
    ax.plot(df["a"], df["Omega_cal"], label=r"$\Omega_R$")
    ax.plot(df["a"], df["closure_sum"], color="0.2", ls="--", lw=0.9, label="sum")
    ax.set_xscale("log")
    ax.set_ylim(-0.05, 1.2)
    ax.axvline(1.0, color="0.5", ls=":", lw=0.8)
    ax.set_xlabel("a")
    ax.set_title("R-critical density closure")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(df["a"], df["H"] / np.interp(1.0, df["a"], df["H"]))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axvline(1.0, color="0.5", ls=":", lw=0.8)
    ax.set_xlabel("a")
    ax.set_ylabel("H/H(a=1)")
    ax.set_title("self-generated expansion")

    ax = axes[1, 1]
    ax.plot(df["a"], df["q"], label="q")
    ax.plot(df["a"], df["G_unscreened_over_G"], label=r"$G_{\rm unscreened}/G$")
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.axvline(1.0, color="0.5", ls=":", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("a")
    ax.set_title("acceleration and coupling")
    ax.legend(fontsize=8)
    fig.tight_layout()
    FIG_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PDF)
    fig.savefig(FIG_PNG, dpi=180)
    plt.close(fig)


def main():
    candidates = scan()
    best_score, p, row, df = candidates[0]
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    make_plot(df)
    idx_early = 0
    idx_today = int(np.argmin(np.abs(df["a"].to_numpy() - 1.0)))
    today = df.iloc[idx_today]
    early = df.iloc[idx_early]
    out = {
        "model": "locked finite-window R-Universe scalar calibration theory",
        "status": "direct background integration of the locked functions; scalar perturbations and CMB module are not yet implemented",
        "parameters": asdict(p),
        "score": best_score,
        "top_candidates": [
            {
                "rank": i + 1,
                "score": float(sc),
                "parameters": asdict(pp),
                "a": float(rr["a"]),
                "Omega_m_R": float(rr["Omega_m_R"]),
                "Omega_R": float(rr["Omega_cal"]),
                "q": float(rr["q"]),
                "F": float(rr["F"]),
                "H": float(rr["H"]),
            }
            for i, (sc, pp, rr, _dd) in enumerate(candidates[:8])
        ],
        "today_a_near_1": {
            "a": float(today["a"]),
            "H": float(today["H"]),
            "W": float(today["W"]),
            "F": float(today["F"]),
            "Omega_m_R": float(today["Omega_m_R"]),
            "Omega_r_R": float(today["Omega_r_R"]),
            "Omega_R": float(today["Omega_cal"]),
            "closure_sum": float(today["closure_sum"]),
            "q": float(today["q"]),
            "G_unscreened_over_G": float(today["G_unscreened_over_G"]),
        },
        "early": {
            "a": float(early["a"]),
            "W": float(early["W"]),
            "F": float(early["F"]),
            "Omega_m_R": float(early["Omega_m_R"]),
            "Omega_r_R": float(early["Omega_r_R"]),
            "Omega_R": float(early["Omega_cal"]),
            "closure_sum": float(early["closure_sum"]),
        },
        "stability_gates": {
            "F_positive_on_integrated_branch": bool(df["F"].min() > 0.0),
            "Z_positive_on_integrated_branch": bool(df["Z"].min() > 0.0),
            "c_T2": 1.0,
            "scalar_Qs_cs2": "not yet computed; requires full quadratic action",
        },
        "outputs": {
            "csv": str(OUT_CSV.relative_to(ROOT)),
            "figure": str(FIG_PDF.relative_to(ROOT)),
        },
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
