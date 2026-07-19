#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from r_universe_core import (
    CC_Z,
    DATA,
    RD_FIXED,
    SN_MB,
    chi2_bao,
    chi2_cc,
    chi2_sn_full,
    coordinate_refine,
    load_bao,
)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "extended_results.json"


def total(x, model, rows, icov):
    if model == "lcdm":
        h0, om = x
        theta, nu = 2.0, 0.0
    elif model == "r1":
        h0, om, theta = x
        nu = 0.0
    else:
        h0, om, theta, nu = x
    if not (50 < h0 < 90 and 0.05 < om < 0.65 and 0.5 < theta < 4.0 and -1.5 < nu < 1.5):
        return 1e99
    try:
        return chi2_sn_full(h0, om, theta, nu) + chi2_bao(h0, om, theta, nu, rows, icov) + chi2_cc(h0, om, theta, nu)
    except Exception:
        return 1e99


def fit_model(model, rows, cov):
    icov = np.linalg.inv(cov)
    fun = lambda v: total(v, model, rows, icov)
    if model == "lcdm":
        starts = [[68.4, 0.31], [67.6, 0.33], [69.0, 0.29]]
        bounds = [(50, 90), (0.05, 0.65)]
        steps = [0.8, 0.012]
    elif model == "r1":
        starts = [[68.0, 0.293, 1.62], [67.5, 0.289, 1.43], [68.5, 0.31, 2.0]]
        bounds = [(50, 90), (0.05, 0.65), (0.5, 4.0)]
        steps = [0.8, 0.012, 0.10]
    else:
        starts = [
            [68.0, 0.293, 1.62, 0.0],
            [68.0, 0.293, 1.40, 0.2],
            [68.0, 0.293, 1.80, -0.2],
            [68.5, 0.31, 2.0, 0.0],
        ]
        bounds = [(50, 90), (0.05, 0.65), (0.5, 4.0), (-1.5, 1.5)]
        steps = [0.8, 0.012, 0.12, 0.08]
    best = (1e99, None)
    for st in starts:
        val, x = coordinate_refine(st, steps, bounds, fun)
        if val < best[0]:
            best = (val, x)
    return best


def summarize(label, rows, cov):
    ndata = len(SN_MB) + len(rows) + len(CC_Z)
    out = {"label": label, "ndata": ndata, "rd_fixed_mpc": RD_FIXED, "models": {}}
    for model, k in [("lcdm", 2), ("r1", 3), ("r2", 4)]:
        chi, x = fit_model(model, rows, cov)
        if model == "lcdm":
            h0, om = x
            theta, nu = 2.0, 0.0
        elif model == "r1":
            h0, om, theta = x
            nu = 0.0
        else:
            h0, om, theta, nu = x
        icov = np.linalg.inv(cov)
        out["models"][model] = {
            "k": k,
            "H0": float(h0),
            "Omega_m0": float(om),
            "theta": float(theta),
            "nu": float(nu),
            "chi2": float(chi),
            "AIC": float(chi + 2 * k),
            "BIC": float(chi + k * math.log(ndata)),
            "chi2_SN_full": chi2_sn_full(h0, om, theta, nu),
            "chi2_BAO": chi2_bao(h0, om, theta, nu, rows, icov),
            "chi2_CC": chi2_cc(h0, om, theta, nu),
        }
    ref = out["models"]["lcdm"]
    out["delta_vs_lcdm"] = {}
    for model in ["r1", "r2"]:
        m = out["models"][model]
        out["delta_vs_lcdm"][model] = {
            "delta_chi2": m["chi2"] - ref["chi2"],
            "delta_AIC": m["AIC"] - ref["AIC"],
            "delta_BIC": m["BIC"] - ref["BIC"],
        }
    return out


def main():
    dr1_rows, dr1_cov = load_bao(DATA / "desi_dr1/desi_2024_bao_mean.txt", DATA / "desi_dr1/desi_2024_bao_cov.txt")
    dr2_rows, dr2_cov = load_bao(
        DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",
        DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt",
    )
    results = [
        summarize("Pantheon+full + DESI_DR1_BAO + CC", dr1_rows, dr1_cov),
        summarize("Pantheon+full + DESI_DR2_BAO + CC", dr2_rows, dr2_cov),
    ]
    OUT.write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
