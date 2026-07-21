#!/usr/bin/env python3
"""Scan one-parameter running R branches against standard competitors.

This is an adversarial model-building audit.  Each candidate has the same
number of fitted late-time parameters as R1 and wCDM: H0, Omega_m0 and one
shape parameter.  The scan is meant to identify whether a simple relational
running law can improve the late-time likelihood without adding CPL-like
freedom.
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
    load_bao,
)

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "code" / "r_running_scan_results.json"
OUT_CSV = ROOT / "tables" / "r_running_scan.csv"


def theta_of_z(z, alpha, law):
    z = np.asarray(z, dtype=float)
    a = 1.0 / (1.0 + z)
    if law == "constant":
        return 2.0 - alpha + 0.0 * z
    if law == "late_a":
        return 2.0 - alpha * a
    if law == "early_z_over_1pz":
        return 2.0 - alpha * z / (1.0 + z)
    if law == "balanced_1_minus_a2":
        return 2.0 - alpha * (1.0 - a * a)
    if law == "log_running":
        return 2.0 - alpha * np.log1p(z) / np.log(4.0)
    raise ValueError(law)


def e_running(z, om, alpha, law):
    z = np.asarray(z, dtype=float)
    ode = 1.0 - om
    if not (0.0 < om < 1.0 and ode > 0.0 and -4.0 < alpha < 4.0):
        return np.full_like(z, np.nan, dtype=float)
    theta = theta_of_z(z, alpha, law)
    if np.any(theta <= 0.05) or np.any(theta > 6.0):
        return np.full_like(z, np.nan, dtype=float)
    m = om * (1.0 + z) ** 3
    amp = ode
    e = np.sqrt(m + amp)

    def f(x):
        return x * x - amp * np.power(x, 2.0 - theta) - m

    for _ in range(35):
        fp = 2.0 * e - amp * (2.0 - theta) * np.power(e, 1.0 - theta)
        step = f(e) / fp
        en = e - step
        bad = (~np.isfinite(en)) | (en <= 0.0) | (~np.isfinite(step))
        e = np.where(bad, 0.5 * (e + np.sqrt(m + amp)), en)
        if np.nanmax(np.abs(step / e)) < 1e-11:
            break
    scaled = np.nanmax(np.abs(f(e)) / (1.0 + m + e * e))
    if not (np.all(np.isfinite(e)) and np.all(e > 0.0) and scaled < 1e-8):
        lo = np.full_like(z, 1e-12)
        hi = np.sqrt(m + amp) + 8.0
        for _ in range(80):
            bad = f(hi) <= 0.0
            if not np.any(bad):
                break
            hi[bad] *= 2.0
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            go_hi = f(mid) > 0.0
            hi = np.where(go_hi, mid, hi)
            lo = np.where(go_hi, lo, mid)
        e = 0.5 * (lo + hi)
    return e


def background_running(h0, om, alpha, law, zpoints):
    zgrid = np.unique(np.concatenate(([0.0], np.asarray(zpoints, dtype=float))))
    e = e_running(zgrid, om, alpha, law)
    if np.any(~np.isfinite(e)) or np.any(e <= 0.0):
        raise ValueError("bad running-R background")
    inv = 1.0 / e
    dz = np.diff(zgrid)
    integ = np.zeros_like(zgrid)
    integ[1:] = np.cumsum(0.5 * dz * (inv[1:] + inv[:-1]))
    dc = (C / h0) * integ
    return zgrid, e, dc


def chi2_sn_running(h0, om, alpha, law):
    zgrid, _, dc = background_running(h0, om, alpha, law, SN_ZHD)
    integral = np.interp(SN_ZHD, zgrid, dc) / (C / h0)
    dl = (C / h0) * (1.0 + SN_ZHEL) * integral
    if np.any(dl <= 0.0):
        return 1e100
    mu = 5.0 * np.log10(dl) + 25.0
    r = SN_MB - mu
    icr = ICOV_SN @ r
    return float(r @ icr - (ONES_SN @ icr) ** 2 / DEN_SN)


def chi2_bao_running(h0, om, alpha, law, rows, icov):
    zvals = np.array([r[0] for r in rows])
    zgrid, egrid, dcgrid = background_running(h0, om, alpha, law, zvals)
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


def chi2_cc_running(h0, om, alpha, law):
    e = e_running(CC_Z, om, alpha, law)
    return float(np.sum(((CC_H - h0 * e) / CC_SIG) ** 2))


def total_running(x, law, rows, icov):
    h0, om, alpha = x
    if not (50.0 < h0 < 90.0 and 0.05 < om < 0.65 and -3.5 < alpha < 3.5):
        return 1e99
    try:
        return (
            chi2_sn_running(h0, om, alpha, law)
            + chi2_bao_running(h0, om, alpha, law, rows, icov)
            + chi2_cc_running(h0, om, alpha, law)
        )
    except Exception:
        return 1e99


def fit_law(law, rows, cov):
    icov = np.linalg.inv(cov)
    fun = lambda x: total_running(x, law, rows, icov)
    starts = [
        [68.0, 0.293, 0.4],
        [68.6, 0.306, 0.0],
        [67.8, 0.285, 0.8],
        [69.0, 0.315, -0.4],
        [67.5, 0.275, 1.2],
    ]
    best = (1e99, None)
    for st in starts:
        val, x = coordinate_refine(st, [0.8, 0.012, 0.12], [(50, 90), (0.05, 0.65), (-3.5, 3.5)], fun, rounds=10)
        if val < best[0]:
            best = (val, x)
    return best


def main():
    rows, cov = load_bao(DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt", DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt")
    ndata = len(SN_MB) + len(rows) + len(CC_Z)
    lcdm_chi2 = chi2_sn_full(68.63125, 0.3059375, 2.0, 0.0) + chi2_bao(
        68.63125, 0.3059375, 2.0, 0.0, rows, np.linalg.inv(cov)
    ) + chi2_cc(68.63125, 0.3059375, 2.0, 0.0)
    laws = ["constant", "late_a", "early_z_over_1pz", "balanced_1_minus_a2", "log_running"]
    records = []
    for law in laws:
        chi2, x = fit_law(law, rows, cov)
        h0, om, alpha = x
        k = 3
        records.append(
            {
                "model": f"Rrun_{law}",
                "law": law,
                "k": k,
                "H0": float(h0),
                "Omega_m0": float(om),
                "alpha": float(alpha),
                "chi2": float(chi2),
                "AIC": float(chi2 + 2 * k),
                "BIC": float(chi2 + k * math.log(ndata)),
                "delta_chi2_vs_lcdm": float(chi2 - lcdm_chi2),
                "delta_AIC_vs_lcdm": float((chi2 + 2 * k) - (lcdm_chi2 + 4)),
                "delta_BIC_vs_lcdm": float((chi2 + k * math.log(ndata)) - (lcdm_chi2 + 2 * math.log(ndata))),
            }
        )
    records.sort(key=lambda row: row["AIC"])
    out = {
        "dataset": "Pantheon+ full covariance + DESI DR2 BAO + cosmic chronometers",
        "ndata": ndata,
        "lcdm_reference_chi2": float(lcdm_chi2),
        "records": records,
        "best_AIC": records[0],
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))
    pd.DataFrame(records).to_csv(OUT_CSV, index=False)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
