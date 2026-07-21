#!/usr/bin/env python3
"""Fit standard late-time dark-energy competitors on the same public likelihood.

The purpose is not to make R1 look good. It is to test whether the observed
late-time improvement is specific to the R-Universe infrared branch or merely
the generic effect of adding dark-energy shape freedom.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from r_universe_core import (
    C,
    CC_H,
    CC_SIG,
    CC_Z,
    DATA,
    DEN_SN,
    ICOV_SN,
    ONES_SN,
    RD_FIXED,
    SN_MB,
    SN_ZHD,
    SN_ZHEL,
    chi2_bao,
    chi2_cc,
    chi2_sn_full,
    coordinate_refine,
    e_general,
    load_bao,
)

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "code" / "competitor_model_results.json"
OUT_CSV = ROOT / "tables" / "competitor_model_comparison.csv"


def e_competitor(z, om, model, w0=-1.0, wa=0.0):
    z = np.asarray(z, dtype=float)
    ode = 1.0 - om
    if not (0.0 < om < 1.0 and ode > 0):
        return np.full_like(z, np.nan, dtype=float)
    if model == "lcdm":
        de = ode
    elif model == "wcdm":
        de = ode * (1.0 + z) ** (3.0 * (1.0 + w0))
    elif model == "cpl":
        de = ode * (1.0 + z) ** (3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * z / (1.0 + z))
    else:
        raise ValueError(model)
    e2 = om * (1.0 + z) ** 3 + de
    if np.any(e2 <= 0) or np.any(~np.isfinite(e2)):
        return np.full_like(z, np.nan, dtype=float)
    return np.sqrt(e2)


def background_competitor(h0, om, model, w0, wa, zpoints):
    zgrid = np.unique(np.concatenate(([0.0], np.asarray(zpoints, dtype=float))))
    e = e_competitor(zgrid, om, model, w0, wa)
    if np.any(~np.isfinite(e)) or np.any(e <= 0):
        raise ValueError("bad competitor background")
    inv = 1.0 / e
    dz = np.diff(zgrid)
    integ = np.zeros_like(zgrid)
    integ[1:] = np.cumsum(0.5 * dz * (inv[1:] + inv[:-1]))
    dc = (C / h0) * integ
    return zgrid, e, dc


def chi2_sn_competitor(h0, om, model, w0, wa):
    zgrid, _, dc = background_competitor(h0, om, model, w0, wa, SN_ZHD)
    integral = np.interp(SN_ZHD, zgrid, dc) / (C / h0)
    dl = (C / h0) * (1.0 + SN_ZHEL) * integral
    if np.any(dl <= 0):
        return 1e100
    mu = 5.0 * np.log10(dl) + 25.0
    r = SN_MB - mu
    icr = ICOV_SN @ r
    return float(r @ icr - (ONES_SN @ icr) ** 2 / DEN_SN)


def chi2_bao_competitor(h0, om, model, w0, wa, rows, icov):
    zvals = np.array([r[0] for r in rows])
    zgrid, egrid, dcgrid = background_competitor(h0, om, model, w0, wa, zvals)
    pred = []
    for z, _, typ in rows:
        e = float(np.interp(z, zgrid, egrid))
        dm = float(np.interp(z, zgrid, dcgrid))
        dh = C / (h0 * e)
        if typ == "DM_over_rs":
            pred.append(dm / RD_FIXED)
        elif typ == "DH_over_rs":
            pred.append(dh / RD_FIXED)
        elif typ == "DV_over_rs":
            pred.append((z * dm * dm * dh) ** (1.0 / 3.0) / RD_FIXED)
        else:
            raise ValueError(typ)
    resid = np.array([r[1] for r in rows]) - np.array(pred)
    return float(resid @ icov @ resid)


def chi2_cc_competitor(h0, om, model, w0, wa):
    e = e_competitor(CC_Z, om, model, w0, wa)
    return float(np.sum(((CC_H - h0 * e) / CC_SIG) ** 2))


def total_competitor(x, model, rows, icov):
    if model == "lcdm":
        h0, om = x
        w0, wa = -1.0, 0.0
    elif model == "wcdm":
        h0, om, w0 = x
        wa = 0.0
    elif model == "cpl":
        h0, om, w0, wa = x
    else:
        raise ValueError(model)
    if not (50.0 < h0 < 90.0 and 0.05 < om < 0.65 and -2.5 < w0 < -0.2 and -3.0 < wa < 3.0):
        return 1e99
    try:
        return (
            chi2_sn_competitor(h0, om, model, w0, wa)
            + chi2_bao_competitor(h0, om, model, w0, wa, rows, icov)
            + chi2_cc_competitor(h0, om, model, w0, wa)
        )
    except Exception:
        return 1e99


def fit_competitor(model, rows, cov):
    icov = np.linalg.inv(cov)
    fun = lambda v: total_competitor(v, model, rows, icov)
    if model == "lcdm":
        starts = [[68.6, 0.306], [67.5, 0.32], [70.0, 0.28]]
        bounds = [(50, 90), (0.05, 0.65)]
        steps = [0.8, 0.012]
        k = 2
    elif model == "wcdm":
        starts = [[68.6, 0.306, -1.0], [68.0, 0.29, -0.8], [67.5, 0.32, -1.2], [70.0, 0.25, -0.6]]
        bounds = [(50, 90), (0.05, 0.65), (-2.5, -0.2)]
        steps = [0.8, 0.012, 0.10]
        k = 3
    elif model == "cpl":
        starts = [
            [68.6, 0.306, -1.0, 0.0],
            [68.0, 0.29, -0.8, -0.4],
            [67.5, 0.32, -1.2, 0.5],
            [70.0, 0.25, -0.6, -1.0],
        ]
        bounds = [(50, 90), (0.05, 0.65), (-2.5, -0.2), (-3.0, 3.0)]
        steps = [0.8, 0.012, 0.10, 0.20]
        k = 4
    else:
        raise ValueError(model)
    best = (1e99, None)
    for st in starts:
        val, x = coordinate_refine(st, steps, bounds, fun, rounds=8)
        if val < best[0]:
            best = (val, x)
    return k, best[0], best[1]


def r1_from_existing(rows, cov):
    icov = np.linalg.inv(cov)
    candidates = [[67.975, 0.29275, 1.61484375], [68.0, 0.293, 1.62], [68.6, 0.306, 2.0]]

    def objective(x):
        h0, om, theta = x
        if not (50 < h0 < 90 and 0.05 < om < 0.65 and 0.5 < theta < 4.0):
            return 1e99
        try:
            return chi2_sn_full(h0, om, theta, 0.0) + chi2_bao(h0, om, theta, 0.0, rows, icov) + chi2_cc(h0, om, theta, 0.0)
        except Exception:
            return 1e99

    best = (1e99, None)
    for st in candidates:
        val, x = coordinate_refine(st, [0.8, 0.012, 0.10], [(50, 90), (0.05, 0.65), (0.5, 4.0)], objective, rounds=8)
        if val < best[0]:
            best = (val, x)
    return 3, best[0], best[1]


def summarize():
    rows, cov = load_bao(DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt", DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt")
    ndata = len(SN_MB) + len(rows) + len(CC_Z)
    out = {
        "dataset": "Pantheon+ full covariance + DESI DR2 BAO + cosmic chronometers",
        "ndata": ndata,
        "rd_fixed_mpc": RD_FIXED,
        "models": {},
    }
    for model in ["lcdm", "wcdm", "cpl"]:
        k, chi2, x = fit_competitor(model, rows, cov)
        if model == "lcdm":
            h0, om = x
            pars = {"H0": h0, "Omega_m0": om, "w0": -1.0, "wa": 0.0}
        elif model == "wcdm":
            h0, om, w0 = x
            pars = {"H0": h0, "Omega_m0": om, "w0": w0, "wa": 0.0}
        else:
            h0, om, w0, wa = x
            pars = {"H0": h0, "Omega_m0": om, "w0": w0, "wa": wa}
        out["models"][model] = {
            "k": k,
            **{kk: float(vv) for kk, vv in pars.items()},
            "theta": None,
            "chi2": float(chi2),
            "AIC": float(chi2 + 2 * k),
            "BIC": float(chi2 + k * math.log(ndata)),
        }
    k, chi2, x = r1_from_existing(rows, cov)
    h0, om, theta = x
    out["models"]["r1"] = {
        "k": k,
        "H0": float(h0),
        "Omega_m0": float(om),
        "w0": None,
        "wa": None,
        "theta": float(theta),
        "chi2": float(chi2),
        "AIC": float(chi2 + 2 * k),
        "BIC": float(chi2 + k * math.log(ndata)),
    }
    ref = out["models"]["lcdm"]
    for model, m in out["models"].items():
        m["delta_chi2_vs_lcdm"] = float(m["chi2"] - ref["chi2"])
        m["delta_AIC_vs_lcdm"] = float(m["AIC"] - ref["AIC"])
        m["delta_BIC_vs_lcdm"] = float(m["BIC"] - ref["BIC"])
    best_aic = min(out["models"], key=lambda name: out["models"][name]["AIC"])
    best_bic = min(out["models"], key=lambda name: out["models"][name]["BIC"])
    out["selection"] = {"best_AIC": best_aic, "best_BIC": best_bic}
    return out


def main():
    out = summarize()
    OUT_JSON.write_text(json.dumps(out, indent=2))
    pd.DataFrame([{"model": name, **vals} for name, vals in out["models"].items()]).to_csv(OUT_CSV, index=False)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
