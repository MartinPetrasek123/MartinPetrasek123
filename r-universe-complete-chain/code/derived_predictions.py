#!/usr/bin/env python3
"""Derived R-Universe predictions used by the manuscript.

This script is intentionally self-contained except for the already generated
extended_results.json file. It computes background derivatives, cosmographic
diagnostics, Sandage-Loeb drift, a minimal GR linear-growth solution and an ISW
source proxy for the best Pantheon+ full + DESI DR2 + CC fits.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

C = 299792.458
MPC_KM = 3.0856775814913673e19
YR_S = 365.25 * 24 * 3600
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "code" / "derived_predictions.json"
FIT_JSON = ROOT / "code" / "extended_results.json"


def e_general(z, om, theta=2.0, nu=0.0):
    z = np.asarray(z, dtype=float)
    ode = 1.0 - om
    m = om * (1.0 + z) ** 3
    amp = ode * (1.0 + z) ** nu
    if abs(theta - 2.0) < 1e-12 and abs(nu) < 1e-12:
        return np.sqrt(m + ode)
    e = np.sqrt(m + np.maximum(amp, 1e-30))

    def f(x):
        return x * x - amp * np.power(x, 2.0 - theta) - m

    for _ in range(30):
        fp = 2.0 * e - amp * (2.0 - theta) * np.power(e, 1.0 - theta)
        step = f(e) / fp
        en = e - step
        bad = (~np.isfinite(en)) | (en <= 0) | (~np.isfinite(step))
        e = np.where(bad, 0.7 * e + 0.3 * np.sqrt(m + np.maximum(amp, 1e-30)), en)
        if np.nanmax(np.abs(step / e)) < 1e-11:
            break
    return e


def deriv_z(fun, z, h=1e-4):
    z = float(z)
    if z < h:
        # Five-point forward derivative, needed for stable cosmography at z=0.
        return (-25 * fun(z) + 48 * fun(z + h) - 36 * fun(z + 2 * h) + 16 * fun(z + 3 * h) - 3 * fun(z + 4 * h)) / (12 * h)
    return (fun(z + h) - fun(z - h)) / (2 * h)


def q_of_z(z, om, theta, nu):
    e = lambda x: float(e_general(np.array([x]), om, theta, nu)[0])
    return -1.0 + (1.0 + z) * deriv_z(e, z) / e(z)


def jerk_of_z(z, om, theta, nu):
    q = lambda x: q_of_z(x, om, theta, nu)
    return q(z) + 2.0 * q(z) * q(z) + (1.0 + z) * deriv_z(q, z)


def w_eff(z, om, theta, nu):
    e = lambda x: float(e_general(np.array([x]), om, theta, nu)[0])
    dlnedz = deriv_z(lambda x: math.log(e(x)), z)
    # Effective density is (1-om)(1+z)^nu E^(2-theta).
    dlnrho_dln1pz = nu + (2.0 - theta) * (1.0 + z) * dlnedz
    return -1.0 + dlnrho_dln1pz / 3.0


def growth_solution(om, theta, nu, n=3000):
    xs = np.linspace(math.log(1e-3), 0.0, n + 1)
    dx = xs[1] - xs[0]
    y = np.zeros((n + 1, 2))
    y[0] = [1e-3, 1e-3]

    def omega_m_a(a):
        z = 1.0 / a - 1.0
        e = float(e_general(np.array([z]), om, theta, nu)[0])
        return om * a ** -3 / (e * e)

    def dlnh_dlna(a):
        h = 1e-4
        return (
            math.log(float(e_general(np.array([1.0 / math.exp(math.log(a) + h) - 1.0]), om, theta, nu)[0]))
            - math.log(float(e_general(np.array([1.0 / math.exp(math.log(a) - h) - 1.0]), om, theta, nu)[0]))
        ) / (2 * h)

    def rhs(x, yy):
        a = math.exp(x)
        d, v = yy
        return np.array([v, -(2.0 + dlnh_dlna(a)) * v + 1.5 * omega_m_a(a) * d])

    for i in range(n):
        x = xs[i]
        yy = y[i]
        k1 = rhs(x, yy)
        k2 = rhs(x + 0.5 * dx, yy + 0.5 * dx * k1)
        k3 = rhs(x + 0.5 * dx, yy + 0.5 * dx * k2)
        k4 = rhs(x + dx, yy + dx * k3)
        y[i + 1] = yy + dx * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return xs, y


def growth_at_z(xs, y, z):
    x = math.log(1.0 / (1.0 + z))
    d = np.interp(x, xs, y[:, 0])
    v = np.interp(x, xs, y[:, 1])
    d0 = y[-1, 0]
    return d / d0, v / d


def main():
    fit = json.loads(FIT_JSON.read_text())[1]["models"]
    lcdm = fit["lcdm"]
    r1 = fit["r1"]
    r2 = fit["r2"]
    models = {"lcdm": lcdm, "r1": r1, "r2": r2}

    out = {"source_fit": "Pantheon+full + DESI_DR2_BAO + CC", "models": {}}
    z_grid = np.linspace(0.0, 3.0, 121)
    for name, m in models.items():
        om, th, nu = m["Omega_m0"], m["theta"], m["nu"]
        xs, yy = growth_solution(om, th, nu)
        dz = []
        fz = []
        phiz = []
        for z in z_grid:
            d, f = growth_at_z(xs, yy, float(z))
            e = float(e_general(np.array([z]), om, th, nu)[0])
            dz.append(d)
            fz.append(f)
            phiz.append(d * (1.0 + z))  # proportional to D/a, normalized to today.
        phiz = np.asarray(phiz)
        isw_proxy = np.gradient(phiz, z_grid)
        out["models"][name] = {
            "H0": m["H0"],
            "Omega_m0": om,
            "theta": th,
            "nu": nu,
            "q0": q_of_z(0.0, om, th, nu),
            "j0": jerk_of_z(0.0, om, th, nu),
            "w_eff_0": w_eff(0.0, om, th, nu),
            "w_eff_z1": w_eff(1.0, om, th, nu),
            "H_infty_over_H0_R1_formula": (1.0 - om) ** (1.0 / th) if abs(nu) < 1e-12 else None,
            "sandage_loeb_delta_v_cm_s_yr_z2": 1.0e5 * C * YR_S * (m["H0"] / MPC_KM) * (1.0 - float(e_general(np.array([2.0]), om, th, nu)[0]) / 3.0),
            "growth_D_z0p5": growth_at_z(xs, yy, 0.5)[0],
            "growth_f_z0": growth_at_z(xs, yy, 0.0)[1],
            "growth_f_z1": growth_at_z(xs, yy, 1.0)[1],
            "isw_proxy_z1": float(np.interp(1.0, z_grid, isw_proxy)),
        }

    # Relative R1 diagnostics against LCDM.
    xs_l, yy_l = growth_solution(lcdm["Omega_m0"], lcdm["theta"], lcdm["nu"])
    xs_r, yy_r = growth_solution(r1["Omega_m0"], r1["theta"], r1["nu"])
    rel = []
    for z in np.linspace(0.0, 3.0, 301):
        d_l, f_l = growth_at_z(xs_l, yy_l, float(z))
        d_r, f_r = growth_at_z(xs_r, yy_r, float(z))
        rel.append((float(z), f_r * d_r / (f_l * d_l) - 1.0, d_r / d_l - 1.0))
    zmax, fdmax, _ = min(rel, key=lambda row: row[1])
    out["r1_relative_to_lcdm"] = {
        "max_negative_fD_fraction": fdmax,
        "z_at_max_negative_fD": zmax,
        "fD_suppression_percent": -100.0 * fdmax,
        "D_ratio_z1_minus_one": next(row[2] for row in rel if abs(row[0] - 1.0) < 1e-12),
    }

    RESULTS.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
