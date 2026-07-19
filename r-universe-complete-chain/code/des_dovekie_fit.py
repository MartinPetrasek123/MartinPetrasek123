#!/usr/bin/env python3
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from r_universe_core import (
    DATA,
    RD_FIXED,
    background,
    chi2_bao,
    chi2_cc,
    coordinate_refine,
    load_bao,
)

ROOT = Path(__file__).resolve().parent
DES = DATA / "des_dovekie"
BAO = DATA / "desi_dr2"
OUT = ROOT / "des_dovekie_results.json"


def load_des_sn():
    lines = (DES / "DES-Dovekie_HD.csv").read_text().splitlines()
    names = None
    data_lines = []
    for line in lines:
        if line.startswith("VARNAMES:"):
            names = line.split()[1:]
        elif line.strip() and not line.startswith("#"):
            if not line.startswith("VARNAMES:"):
                data_lines.append(line)
    if names is None:
        raise RuntimeError("DES-Dovekie VARNAMES header not found")
    from io import StringIO

    hd = pd.read_csv(StringIO("\n".join(data_lines)), sep=r"\s+", names=names)
    d = np.load(DES / "STAT+SYS.npz")
    n = int(d[d.files[0]][0])
    inv = np.zeros((n, n))
    inv[np.triu_indices(n)] = d[d.files[1]]
    low = np.tril_indices(n, -1)
    inv[low] = inv.T[low]
    mask = hd["zHD"].to_numpy(float) > 0.0
    idx = np.where(mask)[0]
    return (
        hd.loc[mask, "zHD"].to_numpy(float),
        hd.loc[mask, "zHEL"].to_numpy(float),
        hd.loc[mask, "MU"].to_numpy(float),
        inv[np.ix_(idx, idx)],
    )


ZHD, ZHEL, MU, INV = load_des_sn()
ONES = np.ones(len(MU))
INV1 = INV @ ONES
DEN = float(ONES @ INV1)


def chi2_des_sn(h0, om, theta, nu):
    zgrid, _, dc = background(h0, om, theta, nu, ZHD)
    dci = np.interp(ZHD, zgrid, dc)
    dl = (1.0 + ZHEL) * dci
    if np.any(dl <= 0):
        return 1e100
    mu_model = 5.0 * np.log10(dl) + 25.0
    # Same analytic offset marginalization as the DES-Dovekie likelihood, but
    # without the additive normalization because model comparisons share it.
    r = mu_model - MU
    ir = INV @ r
    return float(r @ ir - (ONES @ ir) ** 2 / DEN)


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
        return chi2_des_sn(h0, om, theta, nu) + chi2_bao(h0, om, theta, nu, rows, icov) + chi2_cc(h0, om, theta, nu)
    except Exception:
        return 1e99


def fit_model(model, rows, cov):
    icov = np.linalg.inv(cov)
    fun = lambda v: total(v, model, rows, icov)
    if model == "lcdm":
        starts = [[68.6, 0.306], [70.0, 0.30], [67.0, 0.33]]
        bounds = [(50, 90), (0.05, 0.65)]
        steps = [0.8, 0.012]
    elif model == "r1":
        starts = [[68.0, 0.293, 1.61], [70.0, 0.28, 1.3], [68.6, 0.306, 2.0]]
        bounds = [(50, 90), (0.05, 0.65), (0.5, 4.0)]
        steps = [0.8, 0.012, 0.12]
    else:
        starts = [[68.0, 0.293, 1.61, 0.0], [68.0, 0.30, 2.4, 0.5], [70.0, 0.28, 1.2, -0.2]]
        bounds = [(50, 90), (0.05, 0.65), (0.5, 4.0), (-1.5, 1.5)]
        steps = [0.8, 0.012, 0.12, 0.08]
    best = (1e99, None)
    for st in starts:
        val, x = coordinate_refine(st, steps, bounds, fun)
        if val < best[0]:
            best = (val, x)
    return best


def main():
    rows, cov = load_bao(BAO / "desi_gaussian_bao_ALL_GCcomb_mean.txt", BAO / "desi_gaussian_bao_ALL_GCcomb_cov.txt")
    ndata = len(MU) + len(rows) + 31
    out = {"label": "DES-Dovekie SN STAT+SYS + DESI_DR2_BAO + CC", "ndata": ndata, "rd_fixed_mpc": RD_FIXED, "models": {}}
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
            "chi2_SN_full": chi2_des_sn(h0, om, theta, nu),
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
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
