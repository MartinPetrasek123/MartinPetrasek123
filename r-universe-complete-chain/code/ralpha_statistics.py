#!/usr/bin/env python3
"""Laplace posterior, evidence proxy and look-elsewhere calibration for R-alpha."""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, ncx2

from competitor_models import fit_competitor, total_competitor
from r_running_scan import fit_law, total_running
from r_universe_core import CC_Z, DATA, SN_MB, load_bao

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "code" / "ralpha_statistics_results.json"
OUT_CSV = ROOT / "tables" / "ralpha_laplace_posterior.csv"
OUT_PRIOR_CSV = ROOT / "tables" / "ralpha_prior_sensitivity.csv"
FIG = ROOT / "figures" / "fig07_ralpha_corner.pdf"
FIG_PNG = ROOT / "figures" / "fig07_ralpha_corner.png"


def hessian(fun, x, steps):
    x = np.asarray(x, dtype=float)
    steps = np.asarray(steps, dtype=float)
    n = len(x)
    h = np.zeros((n, n))
    f0 = fun(x)
    for i in range(n):
        xp = x.copy(); xp[i] += steps[i]
        xm = x.copy(); xm[i] -= steps[i]
        h[i, i] = (fun(xp) - 2.0 * f0 + fun(xm)) / (steps[i] ** 2)
        for j in range(i + 1, n):
            xpp = x.copy(); xpp[i] += steps[i]; xpp[j] += steps[j]
            xpm = x.copy(); xpm[i] += steps[i]; xpm[j] -= steps[j]
            xmp = x.copy(); xmp[i] -= steps[i]; xmp[j] += steps[j]
            xmm = x.copy(); xmm[i] -= steps[i]; xmm[j] -= steps[j]
            h[i, j] = h[j, i] = (fun(xpp) - fun(xpm) - fun(xmp) + fun(xmm)) / (4.0 * steps[i] * steps[j])
    return h


def laplace_logz(chi2_min, cov, prior_widths):
    k = len(prior_widths)
    sign, logdet = np.linalg.slogdet(cov)
    if sign <= 0:
        return float("nan")
    return float(-0.5 * chi2_min + 0.5 * k * math.log(2.0 * math.pi) + 0.5 * logdet - np.sum(np.log(prior_widths)))


def make_corner(samples, labels):
    fig, axes = plt.subplots(3, 3, figsize=(7.2, 7.2))
    for i in range(3):
        for j in range(3):
            ax = axes[i, j]
            if i == j:
                ax.hist(samples[:, i], bins=36, color="0.25", alpha=0.85)
            elif i > j:
                ax.scatter(samples[:, j], samples[:, i], s=3, alpha=0.08, color="tab:blue", rasterized=True)
            else:
                ax.axis("off")
            if i == 2 and j <= i:
                ax.set_xlabel(labels[j])
            if j == 0 and i > 0:
                ax.set_ylabel(labels[i])
    fig.suptitle(r"R$\alpha$ Laplace posterior approximation", y=0.94)
    fig.tight_layout()
    fig.savefig(FIG)
    fig.savefig(FIG_PNG, dpi=180)
    plt.close(fig)


def main():
    rows, cov_bao = load_bao(DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt", DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt")
    icov_bao = np.linalg.inv(cov_bao)
    ndata = len(SN_MB) + len(rows) + len(CC_Z)

    ralpha_chi2, ralpha_x = fit_law("late_a", rows, cov_bao)
    lcdm_k, lcdm_chi2, lcdm_x = fit_competitor("lcdm", rows, cov_bao)
    wcdm_k, wcdm_chi2, wcdm_x = fit_competitor("wcdm", rows, cov_bao)

    ralpha_fun = lambda x: total_running(x, "late_a", rows, icov_bao)
    lcdm_fun = lambda x: total_competitor(x, "lcdm", rows, icov_bao)
    wcdm_fun = lambda x: total_competitor(x, "wcdm", rows, icov_bao)

    h_ra = hessian(ralpha_fun, ralpha_x, [0.08, 0.0015, 0.015])
    h_l = hessian(lcdm_fun, lcdm_x, [0.08, 0.0015])
    h_w = hessian(wcdm_fun, wcdm_x, [0.08, 0.0015, 0.015])
    cov_ra = 2.0 * np.linalg.inv(h_ra)
    cov_l = 2.0 * np.linalg.inv(h_l)
    cov_w = 2.0 * np.linalg.inv(h_w)
    sig_ra = np.sqrt(np.diag(cov_ra))
    corr_ra = cov_ra / np.outer(sig_ra, sig_ra)

    # Conservative broad priors matching the optimizer bounds used in the scans.
    priors_broad = {"H0": 40.0, "Omega_m0": 0.60, "alpha": 7.0, "w0": 2.3}
    priors_mid = {"H0": 20.0, "Omega_m0": 0.30, "alpha": 3.0, "w0": 1.2}
    priors_tight = {"H0": 10.0, "Omega_m0": 0.18, "alpha": 1.4, "w0": 0.7}
    prior_sets = {"broad": priors_broad, "mid": priors_mid, "tight": priors_tight}
    prior_records = []
    for name, prior in prior_sets.items():
        logz_l = laplace_logz(lcdm_chi2, cov_l, [prior["H0"], prior["Omega_m0"]])
        logz_w = laplace_logz(wcdm_chi2, cov_w, [prior["H0"], prior["Omega_m0"], prior["w0"]])
        logz_ra = laplace_logz(ralpha_chi2, cov_ra, [prior["H0"], prior["Omega_m0"], prior["alpha"]])
        prior_records.append(
            {
                "prior_set": name,
                "lnZ_lcdm": logz_l,
                "lnZ_wcdm": logz_w,
                "lnZ_ralpha": logz_ra,
                "lnB_ralpha_lcdm": logz_ra - logz_l,
                "lnB_ralpha_wcdm": logz_ra - logz_w,
            }
        )
    pd.DataFrame(prior_records).to_csv(OUT_PRIOR_CSV, index=False)

    rng = np.random.default_rng(20260721)
    delta = float(lcdm_chi2 - ralpha_chi2)
    nsim = 1_000_000
    null_single = rng.chisquare(df=1, size=nsim)
    null_five = np.max(rng.chisquare(df=1, size=(nsim, 5)), axis=1)
    power = rng.noncentral_chisquare(df=1, nonc=delta, size=nsim)
    calib = {
        "observed_delta_chi2_lcdm_minus_ralpha": delta,
        "single_branch_p_chi2_ge_observed": float(np.mean(null_single >= delta)),
        "five_branch_look_elsewhere_p": float(np.mean(null_five >= delta)),
        "analytic_single_branch_p": float(chi2.sf(delta, 1)),
        "analytic_five_branch_p": float(1.0 - chi2.cdf(delta, 1) ** 5),
        "ralpha_recovery_power_at_same_signal_delta_chi2_gt_2": float(np.mean(power >= 2.0)),
        "ralpha_recovery_power_at_same_signal_delta_chi2_gt_observed": float(np.mean(power >= delta)),
        "simulation_count": nsim,
    }

    rng2 = np.random.default_rng(42)
    samples = rng2.multivariate_normal(ralpha_x, cov_ra, size=20000)
    make_corner(samples, [r"$H_0$", r"$\Omega_m$", r"$\alpha$"])

    params = ["H0", "Omega_m0", "alpha"]
    posterior_rows = []
    for i, name in enumerate(params):
        posterior_rows.append(
            {
                "parameter": name,
                "best": float(ralpha_x[i]),
                "sigma_laplace": float(sig_ra[i]),
                "lo_68": float(ralpha_x[i] - sig_ra[i]),
                "hi_68": float(ralpha_x[i] + sig_ra[i]),
                "lo_95": float(ralpha_x[i] - 1.96 * sig_ra[i]),
                "hi_95": float(ralpha_x[i] + 1.96 * sig_ra[i]),
            }
        )
    pd.DataFrame(posterior_rows).to_csv(OUT_CSV, index=False)

    out = {
        "dataset": "Pantheon+ full covariance + DESI DR2 BAO + cosmic chronometers",
        "ndata": ndata,
        "method": "finite-difference Hessian and Laplace approximation around deterministic best fits",
        "ralpha": {
            "best_fit": {"H0": float(ralpha_x[0]), "Omega_m0": float(ralpha_x[1]), "alpha": float(ralpha_x[2])},
            "chi2": float(ralpha_chi2),
            "AIC": float(ralpha_chi2 + 6.0),
            "BIC": float(ralpha_chi2 + 3.0 * math.log(ndata)),
            "covariance": cov_ra.tolist(),
            "correlation": corr_ra.tolist(),
            "sigma": {p: float(s) for p, s in zip(params, sig_ra)},
        },
        "lcdm": {"best_fit": {"H0": float(lcdm_x[0]), "Omega_m0": float(lcdm_x[1])}, "chi2": float(lcdm_chi2)},
        "wcdm": {"best_fit": {"H0": float(wcdm_x[0]), "Omega_m0": float(wcdm_x[1]), "w0": float(wcdm_x[2])}, "chi2": float(wcdm_chi2)},
        "laplace_evidence_prior_sensitivity": prior_records,
        "mock_delta_chi2_calibration": calib,
        "outputs": {
            "posterior_csv": "tables/ralpha_laplace_posterior.csv",
            "prior_sensitivity_csv": "tables/ralpha_prior_sensitivity.csv",
            "corner_plot": "figures/fig07_ralpha_corner.pdf",
        },
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
