#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from r_universe_core import C, e_general
from derived_predictions import jerk_of_z

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "code" / "numerical_validation.json"


def e_bisection(z, om, theta, nu=0.0):
    z = np.asarray(z, dtype=float)
    m = om * (1 + z) ** 3
    amp = (1 - om) * (1 + z) ** nu

    def f(x):
        return x * x - amp * np.power(x, 2.0 - theta) - m

    lo = np.full_like(z, 1e-14)
    hi = np.sqrt(m + np.maximum(amp, 1e-30)) + 8.0
    for _ in range(60):
        bad = f(hi) <= 0
        if not np.any(bad):
            break
        hi[bad] *= 2
    for _ in range(140):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        hi = np.where(fm > 0, mid, hi)
        lo = np.where(fm > 0, lo, mid)
    return 0.5 * (lo + hi)


def simpson_distance(z, h0, om, theta, nu=0.0, n=4096):
    if n % 2:
        n += 1
    grid = np.linspace(0, z, n + 1)
    inv = 1.0 / e_bisection(grid, om, theta, nu)
    h = z / n
    integral = h / 3 * (inv[0] + inv[-1] + 4 * inv[1:-1:2].sum() + 2 * inv[2:-1:2].sum())
    return C / h0 * integral


def trapezoid_distance(z, h0, om, theta, nu=0.0, n=4096):
    grid = np.linspace(0, z, n + 1)
    inv = 1.0 / e_general(grid, om, theta, nu)
    return C / h0 * np.trapezoid(inv, grid)


def main():
    z = np.linspace(0, 3.0, 601)
    cases = {
        "lcdm_limit": {"om": 0.3059375, "theta": 2.0, "nu": 0.0},
        "main_r1": {"om": 0.29275, "theta": 1.61484375, "nu": 0.0},
        "main_r2": {"om": 0.30446875, "theta": 2.4078125, "nu": 0.508125},
    }
    out = {"cases": {}}
    for name, p in cases.items():
        diag = {}
        e_newton = e_general(z, p["om"], p["theta"], p["nu"], diagnostics=diag)
        e_bis = e_bisection(z, p["om"], p["theta"], p["nu"])
        m = p["om"] * (1 + z) ** 3
        amp = (1 - p["om"]) * (1 + z) ** p["nu"]
        residual = e_newton**2 - amp * e_newton ** (2 - p["theta"]) - m
        out["cases"][name] = {
            "solver_method": diag["method"],
            "newton_iterations": diag["iterations"],
            "max_scaled_implicit_residual": float(np.max(np.abs(residual) / (1 + m + e_newton**2))),
            "max_relative_newton_vs_bisection": float(np.max(np.abs(e_newton / e_bis - 1))),
            "min_E_on_0_3": float(np.min(e_newton)),
            "max_E_on_0_3": float(np.max(e_newton)),
        }

    om = cases["lcdm_limit"]["om"]
    e_lcdm = e_general(z, om, 2.0, 0.0)
    e_exact = np.sqrt(om * (1 + z) ** 3 + 1 - om)
    d_trap = trapezoid_distance(2.0, 68.63125, om, 2.0)
    d_simp = simpson_distance(2.0, 68.63125, om, 2.0)
    out["lcdm_analytic_validation"] = {
        "max_relative_E_error_theta2": float(np.max(np.abs(e_lcdm / e_exact - 1))),
        "j0_lcdm": float(jerk_of_z(0.0, om, 2.0, 0.0)),
        "distance_z2_trapezoid_mpc": float(d_trap),
        "distance_z2_simpson_mpc": float(d_simp),
        "relative_distance_difference": float(d_trap / d_simp - 1),
    }
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
