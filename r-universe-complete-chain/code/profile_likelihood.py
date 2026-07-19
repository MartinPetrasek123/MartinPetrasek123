#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from r_universe_core import DATA, chi2_bao, chi2_cc, chi2_sn_full, coordinate_refine, load_bao

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "code" / "profile_likelihood_results.json"
OUT_CSV = ROOT / "tables" / "profile_likelihood_theta.csv"


def fit_h0_om_for_theta(theta, rows, cov, start=(67.98, 0.293)):
    icov = np.linalg.inv(cov)

    def obj(v):
        h0, om = v
        if not (50 < h0 < 90 and 0.05 < om < 0.65):
            return 1e99
        return chi2_sn_full(h0, om, theta, 0.0) + chi2_bao(h0, om, theta, 0.0, rows, icov) + chi2_cc(h0, om, theta, 0.0)

    val, x = coordinate_refine(start, [0.45, 0.008], [(50, 90), (0.05, 0.65)], obj, rounds=7)
    return val, x


def interval_from_profile(df, threshold):
    sub = df[df["delta_chi2"] <= threshold]
    if sub.empty:
        return None
    return [float(sub["theta"].min()), float(sub["theta"].max())]


def main():
    rows, cov = load_bao(DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt", DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt")
    # Dense around the minimum and including LCDM theta=2 exactly.
    grid = np.unique(np.concatenate([np.linspace(0.9, 2.6, 86), np.array([1.61484375, 2.0])]))
    prof = []
    start = (67.975, 0.29275)
    for theta in grid:
        chi, x = fit_h0_om_for_theta(float(theta), rows, cov, start=start)
        start = tuple(x)
        prof.append({"theta": float(theta), "chi2_profile": float(chi), "H0_profile": float(x[0]), "Omega_m0_profile": float(x[1])})
    df = pd.DataFrame(prof)
    best_idx = int(df["chi2_profile"].idxmin())
    best = df.loc[best_idx]
    df["delta_chi2"] = df["chi2_profile"] - float(best["chi2_profile"])
    df.to_csv(OUT_CSV, index=False)
    lcdm_row = df.iloc[(df["theta"] - 2.0).abs().argmin()]
    out = {
        "dataset": "Pantheon+full + DESI_DR2_BAO + CC",
        "method": "profile likelihood over theta; H0 and Omega_m0 re-optimized at each fixed theta",
        "theta_grid_min": float(grid.min()),
        "theta_grid_max": float(grid.max()),
        "n_grid": int(len(grid)),
        "best": {k: float(best[k]) for k in best.index},
        "theta_interval_delta_chi2_le_1": interval_from_profile(df, 1.0),
        "theta_interval_delta_chi2_le_4": interval_from_profile(df, 4.0),
        "lcdm_theta2_profile": {k: float(lcdm_row[k]) for k in lcdm_row.index},
        "likelihood_ratio_sigma_approx_theta2": float(math.sqrt(max(float(lcdm_row["delta_chi2"]), 0.0))),
        "profile_csv": str(OUT_CSV.relative_to(ROOT)),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
