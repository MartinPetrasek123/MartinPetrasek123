#!/usr/bin/env python3
"""Evaluate a late-time likelihood from the exact covariant KGB background.

The calculation combines three independent, explicitly supplied data blocks:
Pantheon+ relative distances with the released full covariance and analytic
intercept marginalization, DESI DR2 Gaussian BAO with its full covariance, and
a diagonal cosmic-chronometer H(z) compilation.  The BAO sound horizon is not
fitted or inserted as a fixed number: it is parsed from the H--EFTCAMB run for
the same KGB action and physical matter inputs.

This is a late-time geometric likelihood.  It is not a CMB posterior, does not
include RSD or weak-lensing data, and does not replace a joint parameter
inference with all Planck nuisance parameters.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import cho_factor, cho_solve

from generate_heftcamb_rph import radiation_h
from ru_kgb import RUKGBParams, solve_E_array


ROOT = Path(__file__).resolve().parents[1]
RDRAG_PATTERN = re.compile(r"r_s\(zdrag\)/Mpc\s*=\s*([-+0-9.eE]+)")
C_KM_S = 299792.458


def parse_args() -> argparse.Namespace:
    defaults = RUKGBParams()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True, help="Directory containing pantheon_plus, desi_dr2, and chronometers.")
    parser.add_argument("--solver-dir", type=Path, required=True, help="Completed H--EFTCAMB output directory for this exact parameter point.")
    parser.add_argument("--omega-m0", type=float, default=defaults.omega_m0)
    parser.add_argument("--omega-r0", type=float, default=defaults.omega_r0)
    parser.add_argument("--alpha", type=float, default=defaults.alpha)
    parser.add_argument("--integration-nodes", type=int, default=16385)
    parser.add_argument("--sn-z-min", type=float, default=0.01)
    parser.add_argument("--output", type=Path, default=ROOT / "generated" / "late_time_likelihood.json")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_file(path: Path, description: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"{description} is not available: {path}")
    return path


def parse_rdrag(solver_dir: Path) -> float:
    log = require_file(solver_dir / "solver.log", "H--EFTCAMB solver log")
    values = RDRAG_PATTERN.findall(log.read_text(encoding="utf-8", errors="replace"))
    if len(values) != 1:
        raise RuntimeError(f"expected exactly one r_s(zdrag) entry in {log}, found {len(values)}")
    rdrag = float(values[0])
    if not math.isfinite(rdrag) or rdrag <= 0.0:
        raise RuntimeError("H--EFTCAMB returned an invalid drag sound horizon")
    return rdrag


def load_pantheon(data_root: Path, z_min: float) -> dict[str, Any]:
    directory = data_root / "pantheon_plus"
    table_path = require_file(directory / "pantheon_plus.dat", "Pantheon+ table")
    covariance_path = require_file(directory / "pantheon_plus_stat_sys.cov", "Pantheon+ covariance")
    table = np.genfromtxt(table_path, names=True, dtype=None, encoding="ascii")
    with covariance_path.open(encoding="ascii") as handle:
        size = int(handle.readline().strip())
    covariance = np.loadtxt(covariance_path, skiprows=1).reshape((size, size))
    if len(table) != size:
        raise RuntimeError("Pantheon+ table and covariance dimensions disagree")
    mask = (table["IS_CALIBRATOR"] == 0) & (table["zHD"] > z_min)
    indices = np.flatnonzero(mask)
    if len(indices) == 0:
        raise RuntimeError("Pantheon+ selection is empty")
    covariance = covariance[np.ix_(indices, indices)]
    factor = cho_factor(covariance, lower=True, check_finite=True)
    ones = np.ones(len(indices))
    cinv_ones = cho_solve(factor, ones, check_finite=True)
    return {
        "z_hd": np.asarray(table["zHD"][indices], dtype=float),
        "z_hel": np.asarray(table["zHEL"][indices], dtype=float),
        "m_b_corr": np.asarray(table["m_b_corr"][indices], dtype=float),
        "factor": factor,
        "ones": ones,
        "cinv_ones": cinv_ones,
        "ones_cinv_ones": float(ones @ cinv_ones),
        "count": int(len(indices)),
        "table_sha256": sha256(table_path),
        "covariance_sha256": sha256(covariance_path),
    }


def load_bao(data_root: Path) -> dict[str, Any]:
    directory = data_root / "desi_dr2"
    mean_path = require_file(directory / "desi_gaussian_bao_ALL_GCcomb_mean.txt", "DESI DR2 BAO mean vector")
    covariance_path = require_file(directory / "desi_gaussian_bao_ALL_GCcomb_cov.txt", "DESI DR2 BAO covariance")
    rows: list[tuple[float, float, str]] = []
    for line in mean_path.read_text(encoding="ascii").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        z, value, quantity = line.split()
        rows.append((float(z), float(value), quantity))
    covariance = np.loadtxt(covariance_path)
    if covariance.shape != (len(rows), len(rows)):
        raise RuntimeError("DESI DR2 BAO covariance dimensions disagree with the mean vector")
    return {
        "rows": rows,
        "factor": cho_factor(covariance, lower=True, check_finite=True),
        "count": len(rows),
        "mean_sha256": sha256(mean_path),
        "covariance_sha256": sha256(covariance_path),
    }


def load_chronometers(data_root: Path) -> dict[str, Any]:
    path = require_file(data_root / "chronometers" / "cosmic_chronometers.csv", "cosmic-chronometer table")
    rows: list[tuple[float, float, float]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            systematic = float(row["sigma_sys"]) if row["sigma_sys"].strip() else 0.0
            sigma = math.hypot(float(row["sigma_stat"]), systematic)
            if sigma <= 0.0:
                raise RuntimeError("cosmic-chronometer table contains a non-positive uncertainty")
            rows.append((float(row["z"]), float(row["H"]), sigma))
    if not rows:
        raise RuntimeError("cosmic-chronometer table is empty")
    return {
        "z": np.array([row[0] for row in rows]),
        "h": np.array([row[1] for row in rows]),
        "sigma": np.array([row[2] for row in rows]),
        "count": len(rows),
        "sha256": sha256(path),
    }


def e_of_z(z: np.ndarray, params: RUKGBParams) -> np.ndarray:
    return solve_E_array(1.0 / (1.0 + z), params)


def distance_grid(params: RUKGBParams, z_max: float, nodes: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if nodes < 1025:
        raise ValueError("at least 1025 integration nodes are required")
    z = np.linspace(0.0, z_max, nodes)
    e = e_of_z(z, params)
    if not np.all(np.isfinite(e)) or np.any(e <= 0.0):
        raise RuntimeError("exact KGB background returned an invalid Hubble function")
    dc_over_hubble_distance = cumulative_trapezoid(1.0 / e, z, initial=0.0)
    return z, e, dc_over_hubble_distance


def evaluate_components(
    params: RUKGBParams,
    rdrag_mpc: float,
    pantheon: dict[str, Any],
    bao: dict[str, Any],
    chronometers: dict[str, Any],
    integration_nodes: int,
) -> dict[str, float]:
    h0 = 100.0 * radiation_h(params)
    z_bao = np.array([row[0] for row in bao["rows"]])
    z_all = np.concatenate((pantheon["z_hd"], z_bao, chronometers["z"]))
    z_grid, e_grid, dc_grid = distance_grid(params, float(np.max(z_all)), integration_nodes)

    dc_sn = np.interp(pantheon["z_hd"], z_grid, dc_grid)
    # The arbitrary Hubble-distance factor is absorbed by the analytically
    # marginalized Pantheon+ intercept, so this retains only physical shape.
    dl_sn = (1.0 + pantheon["z_hel"]) * dc_sn
    if np.any(dl_sn <= 0.0):
        raise RuntimeError("non-positive Pantheon+ luminosity distance")
    sn_residual = pantheon["m_b_corr"] - 5.0 * np.log10(dl_sn)
    cinv_residual = cho_solve(pantheon["factor"], sn_residual, check_finite=True)
    chi2_sn = float(
        sn_residual @ cinv_residual
        - (pantheon["ones"] @ cinv_residual) ** 2 / pantheon["ones_cinv_ones"]
    )

    predictions: list[float] = []
    measurements: list[float] = []
    for redshift, value, quantity in bao["rows"]:
        e = float(np.interp(redshift, z_grid, e_grid))
        dc_mpc = C_KM_S * float(np.interp(redshift, z_grid, dc_grid)) / h0
        dh_mpc = C_KM_S / (h0 * e)
        if quantity == "DM_over_rs":
            prediction = dc_mpc / rdrag_mpc
        elif quantity == "DH_over_rs":
            prediction = dh_mpc / rdrag_mpc
        elif quantity == "DV_over_rs":
            prediction = (redshift * dc_mpc * dc_mpc * dh_mpc) ** (1.0 / 3.0) / rdrag_mpc
        else:
            raise RuntimeError(f"unsupported DESI DR2 BAO quantity: {quantity}")
        predictions.append(prediction)
        measurements.append(value)
    bao_residual = np.asarray(measurements) - np.asarray(predictions)
    chi2_bao = float(bao_residual @ cho_solve(bao["factor"], bao_residual, check_finite=True))

    e_cc = np.interp(chronometers["z"], z_grid, e_grid)
    chi2_cc = float(np.sum(((chronometers["h"] - h0 * e_cc) / chronometers["sigma"]) ** 2))
    return {
        "pantheon_plus_relative": chi2_sn,
        "desi_dr2_bao": chi2_bao,
        "cosmic_chronometers": chi2_cc,
        "total": chi2_sn + chi2_bao + chi2_cc,
        "H0_km_s_Mpc": h0,
        "r_drag_Mpc": rdrag_mpc,
    }


def main() -> None:
    args = parse_args()
    params = RUKGBParams(omega_m0=args.omega_m0, omega_r0=args.omega_r0, alpha=args.alpha)
    params.validate()
    rdrag_mpc = parse_rdrag(args.solver_dir)
    pantheon = load_pantheon(args.data_root, args.sn_z_min)
    bao = load_bao(args.data_root)
    chronometers = load_chronometers(args.data_root)
    components = evaluate_components(params, rdrag_mpc, pantheon, bao, chronometers, args.integration_nodes)
    output: dict[str, Any] = {
        "model": "R-Universe covariant KGB exact background",
        "scope": (
            "Late-time geometric likelihood with Pantheon+ relative distances, DESI DR2 BAO, and cosmic chronometers. "
            "The sound horizon is read from the native H--EFTCAMB run of this exact point. "
            "This is not a CMB posterior, RSD/weak-lensing likelihood, or full multi-probe model comparison."
        ),
        "parameters": {
            "omega_m0": params.omega_m0,
            "omega_r0": params.omega_r0,
            "omega_R0": params.omega_R0,
            "alpha": params.alpha,
        },
        "solver": {
            "directory": str(args.solver_dir),
            "solver_log_sha256": sha256(require_file(args.solver_dir / "solver.log", "H--EFTCAMB solver log")),
            "r_drag_source": "H--EFTCAMB solver.log r_s(zdrag)/Mpc",
        },
        "data": {
            "pantheon_plus": {
                "count": pantheon["count"],
                "selection": f"IS_CALIBRATOR=0 and zHD>{args.sn_z_min:g}; additive intercept analytically marginalized",
                "table_sha256": pantheon["table_sha256"],
                "covariance_sha256": pantheon["covariance_sha256"],
            },
            "desi_dr2_bao": {
                "count": bao["count"],
                "mean_sha256": bao["mean_sha256"],
                "covariance_sha256": bao["covariance_sha256"],
            },
            "cosmic_chronometers": {
                "count": chronometers["count"],
                "uncertainty": "quadrature of sigma_stat and supplied sigma_sys; diagonal covariance",
                "sha256": chronometers["sha256"],
            },
        },
        "numerics": {"distance_integration_nodes": args.integration_nodes},
        "chi2": components,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
