#!/usr/bin/env python3
"""Shared R-Universe likelihood and background utilities.

All paths are relative to the released package root. This file intentionally
contains no user-specific absolute paths.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

C = 299792.458
RD_FIXED = 147.09
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw"


def load_pantheon_subset():
    sn_dir = DATA / "pantheon_plus"
    df = pd.read_csv(sn_dir / "pantheon_plus.dat", sep=r"\s+")
    n = int((sn_dir / "pantheon_plus_stat_sys.cov").open().readline().strip())
    cov_vals = np.loadtxt(sn_dir / "pantheon_plus_stat_sys.cov", skiprows=1)
    cov = cov_vals.reshape((n, n))
    mask = (df["IS_CALIBRATOR"].to_numpy() == 0) & (df["zHD"].to_numpy() > 0.01)
    idx = np.where(mask)[0]
    df = df.iloc[idx].reset_index(drop=True)
    cov = cov[np.ix_(idx, idx)]
    return df, cov


def load_bao(path_mean: Path, path_cov: Path):
    rows = []
    for line in Path(path_mean).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        z, val, typ = line.split()
        rows.append((float(z), float(val), typ))
    cov = np.loadtxt(path_cov)
    return rows, cov


def load_cc():
    arr = np.array([
        [0.070, 69.0, 19.6],
        [0.090, 69.0, 12.0],
        [0.120, 68.6, 26.2],
        [0.170, 83.0, 8.0],
        [0.179, 75.0, 4.0],
        [0.199, 75.0, 5.0],
        [0.200, 72.9, 29.6],
        [0.270, 77.0, 14.0],
        [0.280, 88.8, 36.6],
        [0.352, 83.0, 14.0],
        [0.3802, 83.0, 13.5],
        [0.400, 95.0, 17.0],
        [0.4004, 77.0, 10.2],
        [0.4247, 87.1, 11.2],
        [0.4497, 92.8, 12.9],
        [0.470, 89.0, 49.6],
        [0.4783, 80.9, 9.0],
        [0.480, 97.0, 62.0],
        [0.593, 104.0, 13.0],
        [0.680, 92.0, 8.0],
        [0.781, 105.0, 12.0],
        [0.875, 125.0, 17.0],
        [0.880, 90.0, 40.0],
        [0.900, 117.0, 23.0],
        [1.037, 154.0, 20.0],
        [1.300, 168.0, 17.0],
        [1.363, 160.0, 33.6],
        [1.430, 177.0, 18.0],
        [1.530, 140.0, 14.0],
        [1.750, 202.0, 40.0],
        [1.965, 186.5, 50.4],
    ])
    return arr[:, 0], arr[:, 1], arr[:, 2]


PANTHEON, COV_SN = load_pantheon_subset()
SN_ZHD = PANTHEON["zHD"].to_numpy(float)
SN_ZHEL = PANTHEON["zHEL"].to_numpy(float)
SN_MB = PANTHEON["m_b_corr"].to_numpy(float)
ICOV_SN = np.linalg.inv(COV_SN)
ONES_SN = np.ones(len(PANTHEON))
DEN_SN = float(ONES_SN @ ICOV_SN @ ONES_SN)
CC_Z, CC_H, CC_SIG = load_cc()


def e_general(z, om, theta=2.0, nu=0.0, diagnostics: dict | None = None):
    z = np.asarray(z, dtype=float)
    ode = 1.0 - om
    if om <= 0 or ode <= 0 or theta <= 0:
        return np.full_like(z, np.nan, dtype=float)
    m = om * (1.0 + z) ** 3
    amp = ode * (1.0 + z) ** nu
    if abs(theta - 2.0) < 1e-12 and abs(nu) < 1e-12:
        e = np.sqrt(m + ode)
        if diagnostics is not None:
            diagnostics["method"] = "lcdm_analytic"
            diagnostics["iterations"] = 0
            diagnostics["max_residual"] = 0.0
        return e
    e = np.sqrt(m + np.maximum(amp, 1e-30))

    def f(x):
        return x * x - amp * np.power(x, 2.0 - theta) - m

    method = "newton"
    iterations = 0
    for iterations in range(1, 31):
        fp = 2.0 * e - amp * (2.0 - theta) * np.power(e, 1.0 - theta)
        step = f(e) / fp
        en = e - step
        bad = (~np.isfinite(en)) | (en <= 0) | (~np.isfinite(step))
        e = np.where(bad, 0.7 * e + 0.3 * np.sqrt(m + np.maximum(amp, 1e-30)), en)
        if np.nanmax(np.abs(step / e)) < 1e-11:
            break
    scaled = np.nanmax(np.abs(f(e)) / (1.0 + m + e * e))
    if not (np.all(np.isfinite(e)) and np.all(e > 0) and scaled < 1e-9):
        method = "bisection"
        lo = np.full_like(z, 1e-12)
        hi = np.sqrt(m + np.maximum(amp, 1e-30)) + 8.0
        for _ in range(40):
            bad = f(hi) <= 0
            if not np.any(bad):
                break
            hi[bad] *= 2.0
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            fm = f(mid)
            go_hi = fm > 0
            hi = np.where(go_hi, mid, hi)
            lo = np.where(go_hi, lo, mid)
        e = 0.5 * (lo + hi)
        scaled = np.nanmax(np.abs(f(e)) / (1.0 + m + e * e))
    if diagnostics is not None:
        diagnostics["method"] = method
        diagnostics["iterations"] = int(iterations)
        diagnostics["max_residual"] = float(scaled)
    return e


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
    zgrid, _, dc = background(h0, om, theta, nu, SN_ZHD)
    integral = np.interp(SN_ZHD, zgrid, dc) / (C / h0)
    dl = (C / h0) * (1.0 + SN_ZHEL) * integral
    if np.any(dl <= 0):
        return 1e100
    mu = 5.0 * np.log10(dl) + 25.0
    r = SN_MB - mu
    icr = ICOV_SN @ r
    return float(r @ icr - (ONES_SN @ icr) ** 2 / DEN_SN)


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
    e = e_general(CC_Z, om, theta, nu)
    return float(np.sum(((CC_H - h0 * e) / CC_SIG) ** 2))


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
