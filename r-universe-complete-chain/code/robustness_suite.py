#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from r_universe_core import (
    C,
    CC_Z,
    DATA,
    RD_FIXED,
    SN_MB,
    background,
    chi2_bao,
    chi2_cc,
    chi2_sn_full,
    coordinate_refine,
    load_bao,
)

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "code" / "robustness_results.json"
OUT_CSV = ROOT / "tables" / "robustness_suite.csv"


def chi2_bao_with_rd(h0, om, theta, nu, rd, rows, icov):
    zvals = np.array([r[0] for r in rows])
    zgrid, egrid, dcgrid = background(h0, om, theta, nu, zvals)
    pred = []
    for z, _, typ in rows:
        e = float(np.interp(z, zgrid, egrid))
        dm = float(np.interp(z, zgrid, dcgrid))
        dh = C / (h0 * e)
        if typ == "DM_over_rs":
            pred.append(dm / rd)
        elif typ == "DH_over_rs":
            pred.append(dh / rd)
        elif typ == "DV_over_rs":
            pred.append((z * dm * dm * dh) ** (1 / 3) / rd)
        else:
            raise ValueError(typ)
    resid = np.array([r[1] for r in rows]) - np.array(pred)
    return float(resid @ icov @ resid)


def model_params(x, model, free_rd=False):
    if model == "lcdm":
        h0, om = x[:2]
        theta, nu = 2.0, 0.0
        rd = x[2] if free_rd else RD_FIXED
    else:
        h0, om, theta = x[:3]
        nu = 0.0
        rd = x[3] if free_rd else RD_FIXED
    return h0, om, theta, nu, rd


def fit_case(case_name, use_sn, use_bao, use_cc, rows, cov, free_rd=False):
    icov = np.linalg.inv(cov) if use_bao else None
    ndata = (len(SN_MB) if use_sn else 0) + (len(rows) if use_bao else 0) + (len(CC_Z) if use_cc else 0)

    def objective(x, model):
        h0, om, theta, nu, rd = model_params(x, model, free_rd)
        if not (50 < h0 < 90 and 0.05 < om < 0.65 and 0.5 < theta < 4.0 and 130 < rd < 165):
            return 1e99
        total = 0.0
        if use_sn:
            total += chi2_sn_full(h0, om, theta, nu)
        if use_bao:
            if free_rd:
                total += chi2_bao_with_rd(h0, om, theta, nu, rd, rows, icov)
            else:
                total += chi2_bao(h0, om, theta, nu, rows, icov)
        if use_cc:
            total += chi2_cc(h0, om, theta, nu)
        return total

    out = {"case": case_name, "ndata": ndata, "free_rd": free_rd, "models": {}}
    for model in ["lcdm", "r1"]:
        if model == "lcdm" and free_rd:
            starts = [[68.5, 0.306, 147.09], [67.0, 0.32, 147.0], [70.0, 0.28, 145.0]]
            bounds = [(50, 90), (0.05, 0.65), (130, 165)]
            steps = [0.8, 0.012, 0.6]
            k = 3
        elif model == "lcdm":
            starts = [[68.5, 0.306], [67.0, 0.32], [70.0, 0.28]]
            bounds = [(50, 90), (0.05, 0.65)]
            steps = [0.8, 0.012]
            k = 2
        elif free_rd:
            starts = [[68.0, 0.293, 1.61, 147.09], [67.5, 0.30, 2.0, 146.8], [70.0, 0.28, 1.4, 145.0]]
            bounds = [(50, 90), (0.05, 0.65), (0.5, 4.0), (130, 165)]
            steps = [0.8, 0.012, 0.10, 0.6]
            k = 4
        else:
            starts = [[68.0, 0.293, 1.61], [67.5, 0.30, 2.0], [70.0, 0.28, 1.4]]
            bounds = [(50, 90), (0.05, 0.65), (0.5, 4.0)]
            steps = [0.8, 0.012, 0.10]
            k = 3
        best = (1e99, None)
        for st in starts:
            val, x = coordinate_refine(st, steps, bounds, lambda y, m=model: objective(y, m), rounds=6)
            if val < best[0]:
                best = (val, x)
        h0, om, theta, nu, rd = model_params(best[1], model, free_rd)
        out["models"][model] = {
            "k": k,
            "H0": float(h0),
            "Omega_m0": float(om),
            "theta": float(theta),
            "nu": float(nu),
            "rd": float(rd),
            "chi2": float(best[0]),
            "AIC": float(best[0] + 2 * k),
            "BIC": float(best[0] + k * math.log(max(ndata, 2))),
        }
    ref = out["models"]["lcdm"]
    r1 = out["models"]["r1"]
    out["delta_vs_lcdm"] = {
        "delta_chi2": r1["chi2"] - ref["chi2"],
        "delta_AIC": r1["AIC"] - ref["AIC"],
        "delta_BIC": r1["BIC"] - ref["BIC"],
        "sqrt_delta_chi2_sigma_approx": math.sqrt(max(ref["chi2"] - r1["chi2"], 0.0)),
    }
    return out


def flatten_case(case):
    rows = []
    d = case["delta_vs_lcdm"]
    for model, m in case["models"].items():
        rows.append({
            "case": case["case"],
            "model": model,
            "ndata": case["ndata"],
            "free_rd": case["free_rd"],
            **m,
            "delta_chi2_R1_minus_LCDM": d["delta_chi2"] if model == "r1" else 0.0,
            "delta_AIC_R1_minus_LCDM": d["delta_AIC"] if model == "r1" else 0.0,
            "delta_BIC_R1_minus_LCDM": d["delta_BIC"] if model == "r1" else 0.0,
            "sigma_approx": d["sqrt_delta_chi2_sigma_approx"] if model == "r1" else 0.0,
        })
    return rows


def main():
    rows_bao, cov_bao = load_bao(DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt", DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt")
    cases = [
        fit_case("Pantheon_only", True, False, False, rows_bao, cov_bao),
        fit_case("DESI_DR2_BAO_only", False, True, False, rows_bao, cov_bao),
        fit_case("CC_only", False, False, True, rows_bao, cov_bao),
        fit_case("Pantheon_DESI_DR2_BAO", True, True, False, rows_bao, cov_bao),
        fit_case("Pantheon_DESI_DR2_BAO_CC", True, True, True, rows_bao, cov_bao),
        fit_case("Pantheon_DESI_DR2_BAO_CC_free_rd", True, True, True, rows_bao, cov_bao, free_rd=True),
    ]
    # Leave-one-observable BAO jackknife for the main SN+BAO+CC case.
    for i, row in enumerate(rows_bao):
        mask = [j for j in range(len(rows_bao)) if j != i]
        sub_rows = [rows_bao[j] for j in mask]
        sub_cov = cov_bao[np.ix_(mask, mask)]
        label = f"jackknife_without_BAO_{i:02d}_{row[2]}_z{row[0]:.3f}"
        cases.append(fit_case(label, True, True, True, sub_rows, sub_cov))
    OUT_JSON.write_text(json.dumps(cases, indent=2))
    flat = []
    for case in cases:
        flat.extend(flatten_case(case))
    pd.DataFrame(flat).to_csv(OUT_CSV, index=False)
    print(json.dumps(cases, indent=2))


if __name__ == "__main__":
    main()
