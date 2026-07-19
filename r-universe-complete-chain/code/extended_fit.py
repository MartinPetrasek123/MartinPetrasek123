#!/usr/bin/env python3
import json
import math
import sys
from pathlib import Path

import numpy as np

SRC = Path("/Users/mpetr/Documents/Codex/2026-07-09/v/work/r_universe_fit")
sys.path.insert(0, str(SRC))
import fit_r_universe as base  # noqa: E402

C = 299792.458
RD_FIXED = 147.09
ROOT = Path(__file__).resolve().parent
OUT = ROOT / "extended_results.json"


def load_bao(path_mean, path_cov):
    rows = []
    for line in Path(path_mean).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        z, val, typ = line.split()
        rows.append((float(z), float(val), typ))
    cov = np.loadtxt(path_cov)
    return rows, cov


def e_general(z, om, theta=2.0, nu=0.0):
    z = np.asarray(z, dtype=float)
    ode = 1.0 - om
    if om <= 0 or ode <= 0 or theta <= 0:
        return np.full_like(z, np.nan, dtype=float)
    m = om * (1.0 + z) ** 3
    amp = ode * (1.0 + z) ** nu
    if abs(theta - 2.0) < 1e-12 and abs(nu) < 1e-12:
        return np.sqrt(m + ode)
    e = np.sqrt(m + np.maximum(amp, 1e-30))

    def f(x):
        return x * x - amp * np.power(x, 2.0 - theta) - m

    for _ in range(24):
        fp = 2.0 * e - amp * (2.0 - theta) * np.power(e, 1.0 - theta)
        step = f(e) / fp
        en = e - step
        bad = (~np.isfinite(en)) | (en <= 0) | (~np.isfinite(step))
        e = np.where(bad, 0.7 * e + 0.3 * np.sqrt(m + np.maximum(amp, 1e-30)), en)
        if np.nanmax(np.abs(step / e)) < 1e-10:
            break
    if np.all(np.isfinite(e)) and np.all(e > 0) and np.nanmax(np.abs(f(e)) / (1 + m + e * e)) < 1e-8:
        return e

    lo = np.full_like(z, 1e-12)
    hi = np.sqrt(m + np.maximum(amp, 1e-30)) + 8.0
    for _ in range(30):
        bad = f(hi) <= 0
        if not np.any(bad):
            break
        hi[bad] *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        go_hi = fm > 0
        hi = np.where(go_hi, mid, hi)
        lo = np.where(go_hi, lo, mid)
    return 0.5 * (lo + hi)


def background(h0, om, theta, nu, zpoints):
    zgrid = np.unique(np.concatenate(([0.0], np.asarray(zpoints, dtype=float))))
    e = e_general(zgrid, om, theta, nu)
    if np.any(~np.isfinite(e)) or np.any(e <= 0):
        raise ValueError("bad background")
    inv = 1.0 / e
    integ = np.zeros_like(zgrid)
    dz = np.diff(zgrid)
    integ[1:] = np.cumsum(0.5 * dz * (inv[1:] + inv[:-1]))
    dc = (C / h0) * integ
    return zgrid, e, dc


def chi2_sn_full(h0, om, theta, nu):
    zgrid, _, dc = background(h0, om, theta, nu, base.SN_ZHD)
    integral = np.interp(base.SN_ZHD, zgrid, dc) / (C / h0)
    dl = (C / h0) * (1.0 + base.SN_ZHEL) * integral
    if np.any(dl <= 0):
        return 1e100
    mu = 5.0 * np.log10(dl) + 25.0
    r = base.SN_MB - mu
    icr = base.ICOV_SN @ r
    return float(r @ icr - (base.ONES_SN @ icr) ** 2 / base.DEN_SN)


def chi2_bao(h0, om, theta, nu, rows, icov):
    zvals = np.array([r[0] for r in rows])
    zgrid, egrid, dcgrid = background(h0, om, theta, nu, zvals)
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
            pred.append((z * dm * dm * dh) ** (1 / 3) / RD_FIXED)
        else:
            raise ValueError(typ)
    resid = np.array([r[1] for r in rows]) - np.array(pred)
    return float(resid @ icov @ resid)


def chi2_cc(h0, om, theta, nu):
    e = e_general(base.CC_Z, om, theta, nu)
    return float(np.sum(((base.CC_H - h0 * e) / base.CC_SIG) ** 2))


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


def coordinate_refine(x0, steps, bounds, fun, rounds=8):
    x = np.array(x0, dtype=float)
    step = np.array(steps, dtype=float)
    best = fun(x)
    for _ in range(rounds):
        improved = True
        while improved:
            improved = False
            for j in range(len(x)):
                for sgn in (-1, 1):
                    y = x.copy()
                    y[j] += sgn * step[j]
                    if y[j] < bounds[j][0] or y[j] > bounds[j][1]:
                        continue
                    val = fun(y)
                    if val < best:
                        best, x, improved = val, y, True
        step *= 0.5
    return best, x


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
    ndata = len(base.SN_MB) + len(rows) + len(base.CC_Z)
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
    dr1_rows, dr1_cov = load_bao(SRC / "data/desi_2024_bao_mean.txt", SRC / "data/desi_2024_bao_cov.txt")
    dr2_rows, dr2_cov = load_bao(
        ROOT.parent / "public_bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",
        ROOT.parent / "public_bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt",
    )
    results = [summarize("Pantheon+full + DESI_DR1_BAO + CC", dr1_rows, dr1_cov),
               summarize("Pantheon+full + DESI_DR2_BAO + CC", dr2_rows, dr2_cov)]
    OUT.write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
