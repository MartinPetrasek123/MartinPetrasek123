#!/usr/bin/env python3
"""R-alpha minimal relational-kernel predictions.

The script turns the background-only R-alpha fit into a parameter-free
late-time consistency prediction for mu(k,z), eta(k,z), the GW luminosity
ratio Xi(z), and a scale-dependent growth diagnostic.  No new parameter is
fitted here: H0, Omega_m0 and alpha are read from r_running_scan_results.json.

This is a minimal closure, not a completed covariant action.  It is written as
an explicit falsifiable kernel ansatz so every prediction has a code path.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from r_running_scan import e_running

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "code" / "ralpha_master_predictions.json"
OUT_CSV = ROOT / "tables" / "ralpha_master_predictions.csv"
FIG_PDF = ROOT / "figures" / "fig08_ralpha_master_predictions.pdf"
FIG_PNG = ROOT / "figures" / "fig08_ralpha_master_predictions.png"

C = 299792.458
KSTAR = 0.01  # h/Mpc; large-scale clustering/lensing null-test scale
SIGMA8_LCDM_REF = 0.811


def load_best():
    data = json.loads((ROOT / "code" / "r_running_scan_results.json").read_text())
    best = data["best_AIC"]
    return best["H0"], best["Omega_m0"], best["alpha"], data


def e_ra(z, om, alpha):
    return e_running(np.asarray(z, dtype=float), om, alpha, "late_a")


def coupling_beta(z, alpha):
    a = 1.0 / (1.0 + np.asarray(z, dtype=float))
    return alpha * a / (1.0 + alpha * a)


def k_rel_hmpc(z, h0, om, alpha):
    h = h0 / 100.0
    return (h0 * e_ra(z, om, alpha) / C) / h


def window_kernel(k_hmpc, z, h0, om, alpha):
    kr = k_rel_hmpc(z, h0, om, alpha)
    return kr * kr / (k_hmpc * k_hmpc + kr * kr)


def mu_eta(k_hmpc, z, h0, om, alpha):
    beta = coupling_beta(z, alpha)
    win = window_kernel(k_hmpc, z, h0, om, alpha)
    delta = beta * win
    mu = 1.0 + delta
    eta = 1.0 - 0.5 * delta
    sigma = 0.5 * mu * (1.0 + eta)
    return mu, eta, sigma, delta


def xi_of_z(zvals, h0, om, alpha):
    zvals = np.asarray(zvals, dtype=float)
    grid = np.linspace(0.0, max(float(np.max(zvals)), 1e-6), 1200)
    beta = coupling_beta(grid, alpha)
    # Minimal friction closure: d ln M_*^2 / d ln a = beta_R.
    # Xi = exp[1/2 int_0^z beta_R(z') d ln(1+z')].
    integ = np.zeros_like(grid)
    kernel = beta / (1.0 + grid)
    dz = np.diff(grid)
    integ[1:] = np.cumsum(0.5 * dz * (kernel[1:] + kernel[:-1]))
    return np.exp(0.5 * np.interp(zvals, grid, integ))


def growth_solution(h0, om, alpha, k_hmpc=KSTAR):
    def rhs(x, y):
        # x = ln a, z = exp(-x)-1
        a = math.exp(x)
        z = 1.0 / a - 1.0
        eps = 1e-4
        e0 = float(e_ra([z], om, alpha)[0])
        zp = 1.0 / math.exp(x + eps) - 1.0
        zm = 1.0 / math.exp(x - eps) - 1.0
        dlnh = (math.log(float(e_ra([zp], om, alpha)[0])) - math.log(float(e_ra([zm], om, alpha)[0]))) / (2.0 * eps)
        omega_m_a = om * a ** -3 / (e0 * e0)
        mu = float(mu_eta(k_hmpc, np.array([z]), h0, om, alpha)[0][0])
        d, dp = y
        ddp = -(2.0 + dlnh) * dp + 1.5 * mu * omega_m_a * d
        return [dp, ddp]

    x0 = math.log(1.0 / 101.0)
    x1 = 0.0
    y0 = [math.exp(x0), math.exp(x0)]
    sol = solve_ivp(rhs, (x0, x1), y0, rtol=2e-7, atol=1e-9, dense_output=True, max_step=0.02)
    d0 = sol.y[0, -1]

    def eval_z(z):
        x = math.log(1.0 / (1.0 + z))
        d, dp = sol.sol(x)
        d = float(d / d0)
        f = float(dp / sol.sol(x)[0])
        return d, f

    return eval_z


def main():
    h0, om, alpha, scan = load_best()
    zgrid = np.array([0.0, 0.25, 0.5, 1.0, 1.5, 2.0])
    mu, eta, sigma, delta = mu_eta(KSTAR, zgrid, h0, om, alpha)
    xi = xi_of_z(zgrid, h0, om, alpha)
    grow = growth_solution(h0, om, alpha, KSTAR)
    rows = []
    for z, m, et, sg, de, xx in zip(zgrid, mu, eta, sigma, delta, xi):
        d, f = grow(float(z))
        rows.append(
            {
                "z": float(z),
                "k_h_per_Mpc": KSTAR,
                "mu": float(m),
                "eta": float(et),
                "Sigma_lensing": float(sg),
                "Xi_GW_over_EM": float(xx),
                "D_norm": float(d),
                "f_growth": float(f),
                "fsigma8_over_sigma8_ref": float(f * d),
                "delta_kernel": float(de),
                "consistency_eta_minus_one_over_mu_minus_one": float((et - 1.0) / (m - 1.0)) if abs(m - 1.0) > 1e-15 else None,
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(OUT_CSV, index=False)
    dense_z = np.linspace(0.0, 2.0, 240)
    dense_mu, dense_eta, dense_sigma, _ = mu_eta(KSTAR, dense_z, h0, om, alpha)
    dense_xi = xi_of_z(dense_z, h0, om, alpha)
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.4))
    axes[0].plot(dense_z, dense_xi, color="tab:purple")
    axes[0].set_xlabel("z")
    axes[0].set_ylabel(r"$\Xi=d_L^{\rm GW}/d_L^{\rm EM}$")
    axes[0].set_title("standard sirens")
    axes[1].plot(dense_z, dense_mu - 1.0, label=r"$\mu-1$", color="tab:blue")
    axes[1].plot(dense_z, dense_eta - 1.0, label=r"$\eta-1$", color="tab:orange")
    axes[1].set_xlabel("z")
    axes[1].set_ylabel(r"deviation at $k_*$")
    axes[1].set_title("scalar sector")
    axes[1].legend(fontsize=8)
    ratio = (dense_eta - 1.0) / (dense_mu - 1.0)
    axes[2].plot(dense_z, ratio, color="black")
    axes[2].axhline(-0.5, color="0.55", ls="--", lw=1)
    axes[2].set_xlabel("z")
    axes[2].set_ylabel(r"$(\eta-1)/(\mu-1)$")
    axes[2].set_title("master relation")
    fig.tight_layout()
    fig.savefig(FIG_PDF)
    fig.savefig(FIG_PNG, dpi=180)
    plt.close(fig)
    out = {
        "model": "Ralpha minimal relational-kernel closure",
        "source_fit": "r_running_scan_results.json best_AIC",
        "parameters": {"H0": h0, "Omega_m0": om, "alpha": alpha, "k_star_h_per_Mpc": KSTAR},
        "closure": {
            "theta_a": "theta(a)=2-alpha*a",
            "beta_R": "alpha*a/(1+alpha*a)",
            "k_R_h_per_Mpc": "H(z)/(c h)",
            "mu_minus_1": "beta_R*k_R^2/(k^2+k_R^2)",
            "eta_minus_1": "-0.5*(mu-1)",
            "Xi": "exp[0.5*int_0^z beta_R(z') d ln(1+z')]",
            "master_relation_1": "(eta-1)/(mu-1)=-1/2",
            "master_relation_2": "d ln Xi / d ln(1+z)=beta_R/2",
        },
        "predictions": rows,
        "killer_values": {
            "Xi_z1": float(xi_of_z(np.array([1.0]), h0, om, alpha)[0]),
            "mu_kstar_z1": float(mu_eta(KSTAR, np.array([1.0]), h0, om, alpha)[0][0]),
            "eta_kstar_z1": float(mu_eta(KSTAR, np.array([1.0]), h0, om, alpha)[1][0]),
            "fD_kstar_z1": float(grow(1.0)[0] * grow(1.0)[1]),
        },
        "outputs": {
            "prediction_csv": "tables/ralpha_master_predictions.csv",
            "figure": "figures/fig08_ralpha_master_predictions.pdf",
        },
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
