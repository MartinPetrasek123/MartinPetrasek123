#!/usr/bin/env python3
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1]
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

from r_universe_core import DATA, ICOV_SN, RD_FIXED, SN_MB, SN_ZHD, SN_ZHEL, background, e_general, load_bao  # noqa: E402
from derived_predictions import growth_at_z, growth_solution  # noqa: E402

C = 299792.458


def model_mu(zcmb, zhel, h0, om, theta, nu=0.0):
    zg, _, dc = background(h0, om, theta, nu, zcmb)
    dci = np.interp(zcmb, zg, dc)
    dl = (1.0 + zhel) * dci
    return 5.0 * np.log10(dl) + 25.0


def offset_marginalized_residuals(mu_obs, mu_model, invcov):
    ones = np.ones(len(mu_obs))
    r = mu_obs - mu_model
    # Weighted best additive offset for data - model.
    off = float((ones @ invcov @ r) / (ones @ invcov @ ones))
    return r - off, off


def weighted_bins(x, y, nbins=28):
    order = np.argsort(x)
    x = np.asarray(x)[order]
    y = np.asarray(y)[order]
    chunks = np.array_split(np.arange(len(x)), nbins)
    bx, by, be = [], [], []
    for ch in chunks:
        if len(ch) == 0:
            continue
        bx.append(float(np.mean(x[ch])))
        by.append(float(np.mean(y[ch])))
        be.append(float(np.std(y[ch]) / max(math.sqrt(len(ch)), 1)))
    return np.array(bx), np.array(by), np.array(be)


def bao_prediction(rows, h0, om, theta, nu=0.0):
    zvals = np.array([r[0] for r in rows])
    zg, eg, dcg = background(h0, om, theta, nu, zvals)
    pred = []
    for z, _, typ in rows:
        e = float(np.interp(z, zg, eg))
        dm = float(np.interp(z, zg, dcg))
        dh = C / (h0 * e)
        if typ == "DM_over_rs":
            pred.append(dm / RD_FIXED)
        elif typ == "DH_over_rs":
            pred.append(dh / RD_FIXED)
        elif typ == "DV_over_rs":
            pred.append((z * dm * dm * dh) ** (1.0 / 3.0) / RD_FIXED)
    return np.array(pred)


def om_diag(z, om, theta, nu=0.0):
    e = e_general(np.asarray(z), om, theta, nu)
    return (e * e - 1.0) / ((1.0 + np.asarray(z)) ** 3 - 1.0)


def write_summary_tables():
    ext = json.loads((OUT / "code/extended_results.json").read_text())
    des = json.loads((OUT / "code/des_dovekie_results.json").read_text())
    rows = []
    for block in ext + [des]:
        ref = block["models"]["lcdm"]
        for name, m in block["models"].items():
            rows.append({
                "dataset": block["label"],
                "model": name,
                "ndata": block["ndata"],
                "k": m["k"],
                "H0": m["H0"],
                "Omega_m0": m["Omega_m0"],
                "theta": m["theta"],
                "nu": m["nu"],
                "chi2": m["chi2"],
                "AIC": m["AIC"],
                "BIC": m["BIC"],
                "delta_chi2": m["chi2"] - ref["chi2"],
                "delta_AIC": m["AIC"] - ref["AIC"],
                "delta_BIC": m["BIC"] - ref["BIC"],
                "chi2_SN": m["chi2_SN_full"],
                "chi2_BAO": m["chi2_BAO"],
                "chi2_CC": m["chi2_CC"],
            })
    pd.DataFrame(rows).to_csv(TAB / "all_model_fits.csv", index=False)

    # BAO residual table for the main Pantheon+ DR2 fit.
    main = ext[1]["models"]
    rows_bao, cov = load_bao(
        DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",
        DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt",
    )
    sig = np.sqrt(np.diag(cov))
    data = np.array([r[1] for r in rows_bao])
    out = []
    for model in ["lcdm", "r1", "r2"]:
        m = main[model]
        pred = bao_prediction(rows_bao, m["H0"], m["Omega_m0"], m["theta"], m["nu"])
        for (z, val, typ), p, s in zip(rows_bao, pred, sig):
            out.append({"model": model, "z": z, "quantity": typ, "data": val, "prediction": p, "sigma_diag": s, "pull_diag": (val - p) / s})
    pd.DataFrame(out).to_csv(TAB / "desi_dr2_bao_predictions.csv", index=False)


def make_plots():
    ext = json.loads((OUT / "code/extended_results.json").read_text())
    des = json.loads((OUT / "code/des_dovekie_results.json").read_text())
    main = ext[1]["models"]
    lcdm = main["lcdm"]
    r1 = main["r1"]
    r2 = main["r2"]

    z = np.linspace(0.001, 2.5, 500)
    plt.figure(figsize=(7.2, 4.6))
    for label, m, color in [("LCDM", lcdm, "black"), ("R1", r1, "tab:blue"), ("R2", r2, "tab:orange")]:
        plt.plot(z, e_general(z, m["Omega_m0"], m["theta"], m["nu"]), label=label, color=color)
    plt.xlabel("z")
    plt.ylabel(r"$E(z)=H(z)/H_0$")
    plt.title(r"Expansion history: Pantheon+ full + DESI DR2 + CC")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "fig01_expansion_Ez.pdf")
    plt.savefig(FIG / "fig01_expansion_Ez.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7.2, 4.6))
    zz = np.linspace(0.02, 3.0, 500)
    plt.plot(zz, om_diag(zz, lcdm["Omega_m0"], 2.0, 0.0), color="black", label="LCDM")
    plt.plot(zz, om_diag(zz, r1["Omega_m0"], r1["theta"], 0.0), color="tab:blue", label="R1")
    plt.plot(zz, om_diag(zz, r2["Omega_m0"], r2["theta"], r2["nu"]), color="tab:orange", label="R2")
    plt.xlabel("z")
    plt.ylabel("Om(z)")
    plt.title(r"$Om(z)$ diagnostic")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "fig02_om_diagnostic.pdf")
    plt.savefig(FIG / "fig02_om_diagnostic.png", dpi=180)
    plt.close()

    rows_bao, cov = load_bao(
        DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",
        DATA / "desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt",
    )
    sig = np.sqrt(np.diag(cov))
    data = np.array([r[1] for r in rows_bao])
    labels = [f"{r[2].replace('_over_rs','')}\nz={r[0]:.3g}" for r in rows_bao]
    x = np.arange(len(rows_bao))
    plt.figure(figsize=(9.6, 4.8))
    plt.axhline(0, color="0.4", lw=1)
    for model, dx, color, name in [(lcdm, -0.22, "black", "LCDM"), (r1, 0, "tab:blue", "R1"), (r2, 0.22, "tab:orange", "R2")]:
        pred = bao_prediction(rows_bao, model["H0"], model["Omega_m0"], model["theta"], model["nu"])
        pull = (data - pred) / sig
        plt.scatter(x + dx, pull, s=28, label=name, color=color)
    plt.xticks(x, labels, rotation=45, ha="right", fontsize=8)
    plt.ylabel("diagonal BAO pull")
    plt.title("DESI DR2 BAO residual diagnostics")
    plt.legend(ncol=3)
    plt.tight_layout()
    plt.savefig(FIG / "fig03_desi_dr2_bao_pulls.pdf")
    plt.savefig(FIG / "fig03_desi_dr2_bao_pulls.png", dpi=180)
    plt.close()

    mu_l = model_mu(SN_ZHD, SN_ZHEL, lcdm["H0"], lcdm["Omega_m0"], 2.0, 0.0)
    mu_r = model_mu(SN_ZHD, SN_ZHEL, r1["H0"], r1["Omega_m0"], r1["theta"], 0.0)
    res_l, off_l = offset_marginalized_residuals(SN_MB, mu_l, ICOV_SN)
    res_r, off_r = offset_marginalized_residuals(SN_MB, mu_r, ICOV_SN)
    bx, by_l, be_l = weighted_bins(SN_ZHD, res_l, 34)
    _, by_r, be_r = weighted_bins(SN_ZHD, res_r, 34)
    plt.figure(figsize=(7.4, 4.8))
    plt.axhline(0, color="0.5", lw=1)
    plt.errorbar(bx, by_l, yerr=be_l, fmt="o", ms=3, color="black", label="LCDM")
    plt.errorbar(bx, by_r, yerr=be_r, fmt="o", ms=3, color="tab:blue", label="R1")
    plt.xlabel("zHD")
    plt.ylabel("binned SN residual [mag]")
    plt.title("Pantheon+ residuals after intercept profiling")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "fig04_pantheon_residuals.pdf")
    plt.savefig(FIG / "fig04_pantheon_residuals.png", dpi=180)
    plt.close()

    fit_rows = pd.read_csv(TAB / "all_model_fits.csv")
    sel = fit_rows[(fit_rows["model"].isin(["r1", "r2"]))]
    names = []
    daic = []
    dbic = []
    for _, row in sel.iterrows():
        label = row["dataset"].replace("Pantheon+full + ", "").replace("_BAO", "").replace(" SN STAT+SYS", "")
        names.append(f"{label}\n{row['model'].upper()}")
        daic.append(row["delta_AIC"])
        dbic.append(row["delta_BIC"])
    x = np.arange(len(names))
    plt.figure(figsize=(9.2, 5.2))
    plt.axhline(0, color="0.25", lw=1)
    plt.bar(x - 0.18, daic, width=0.36, label="Delta AIC", color="tab:green")
    plt.bar(x + 0.18, dbic, width=0.36, label="Delta BIC", color="tab:red")
    plt.xticks(x, names, rotation=35, ha="right", fontsize=8)
    plt.ylabel(r"R-Universe minus $\Lambda$CDM")
    plt.title(r"Information criteria: negative values favor R-Universe")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "fig05_model_selection.pdf")
    plt.savefig(FIG / "fig05_model_selection.png", dpi=180)
    plt.close()

    xs_l, yy_l = growth_solution(lcdm["Omega_m0"], lcdm["theta"], lcdm["nu"])
    xs_r, yy_r = growth_solution(r1["Omega_m0"], r1["theta"], r1["nu"])
    z_grid = np.linspace(0.0, 2.5, 80)
    lcdm_d = np.array([growth_at_z(xs_l, yy_l, float(z))[0] for z in z_grid])
    r1_d = np.array([growth_at_z(xs_r, yy_r, float(z))[0] for z in z_grid])
    plt.figure(figsize=(7.2, 4.6))
    plt.plot(z_grid, r1_d, "-", label="R1 diagnostic curve")
    plt.plot(z_grid, lcdm_d, "--", label="LCDM diagnostic curve")
    plt.xlabel("z")
    plt.ylabel(r"$D(z)/D(0)$")
    plt.title("Linear growth normalized to today")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "fig06_growth_Dz.pdf")
    plt.savefig(FIG / "fig06_growth_Dz.png", dpi=180)
    plt.close()


def main():
    write_summary_tables()
    make_plots()
    print("Wrote figures to", FIG)
    print("Wrote tables to", TAB)


if __name__ == "__main__":
    main()
