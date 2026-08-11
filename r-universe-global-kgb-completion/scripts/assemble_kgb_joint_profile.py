#!/usr/bin/env python3
"""Combine executed fixed-input Planck points with an exact KGB late-time scan.

Every point must already have been run by H--EFTCAMB and evaluated by the
official Planck 2018 likelihood.  This utility reads those immutable outputs,
parses the corresponding solver-derived drag sound horizon, and adds the
Pantheon+--DESI DR2 BAO--chronometer likelihood from
``evaluate_kgb_late_time.py``.  It never substitutes a distance prior or a
fixed sound horizon.

The resulting surface is a conditional joint profile: the primordial sector,
Planck calibration, and high-ell nuisance treatment remain fixed to the values
used by the Planck grid.  It is not a posterior or Bayesian evidence result.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from evaluate_kgb_late_time import (
    evaluate_components,
    load_bao,
    load_chronometers,
    load_pantheon,
    parse_rdrag,
)
from ru_kgb import RUKGBParams


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--profile-roots", type=Path, nargs="+", required=True)
    parser.add_argument("--integration-nodes", type=int, default=8193)
    parser.add_argument("--sn-z-min", type=float, default=0.01)
    parser.add_argument("--output", type=Path, default=ROOT / "generated" / "joint_planck_late_time_profile.json")
    return parser.parse_args()


def point_tag(params: RUKGBParams) -> str:
    return f"kgb_alpha_{params.alpha:.8f}_omega_m0_{params.omega_m0:.8f}"


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def records_from_root(root: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    summary_path = root / "summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"missing Planck-profile summary: {summary_path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    records = summary.get("grid", {}).get("records")
    if not isinstance(records, list):
        raise RuntimeError(f"missing grid records in {summary_path}")
    return summary, records


def main() -> None:
    args = parse_args()
    pantheon = load_pantheon(args.data_root, args.sn_z_min)
    bao = load_bao(args.data_root)
    chronometers = load_chronometers(args.data_root)
    unique: dict[tuple[float, float, float], tuple[Path, dict[str, Any]]] = {}
    fixed_inputs: list[dict[str, Any]] = []
    for root in args.profile_roots:
        summary, records = records_from_root(root)
        fixed_inputs.append(summary.get("fixed_inputs", {}))
        for record in records:
            if record.get("status") != "completed":
                continue
            key = (float(record["alpha"]), float(record["omega_m0"]), float(record["omega_r0"]))
            unique.setdefault(key, (root, record))

    assembled: list[dict[str, Any]] = []
    for (alpha, omega_m0, omega_r0), (root, record) in sorted(unique.items()):
        params = RUKGBParams(omega_m0=omega_m0, omega_r0=omega_r0, alpha=alpha)
        params.validate()
        solver_dir = root / point_tag(params)
        planck_path = solver_dir / "planck.json"
        try:
            planck = json.loads(planck_path.read_text(encoding="utf-8"))
            if planck.get("status") != "completed":
                raise RuntimeError("stored Planck evaluation is not completed")
            late_time = evaluate_components(
                params,
                parse_rdrag(solver_dir),
                pantheon,
                bao,
                chronometers,
                args.integration_nodes,
            )
            planck_minus_2logl = float(planck["total"]["minus_2_log_likelihood"])
            assembled.append({
                "status": "completed",
                "alpha": alpha,
                "omega_m0": omega_m0,
                "omega_r0": omega_r0,
                "omega_R0": params.omega_R0,
                "planck_minus_2_log_likelihood": planck_minus_2logl,
                "late_time_chi2": late_time,
                "conditional_joint_statistic": planck_minus_2logl + late_time["total"],
                "solver_directory": relative(solver_dir),
                "planck_record": relative(planck_path),
            })
        except Exception as exc:
            assembled.append({
                "status": "failed",
                "alpha": alpha,
                "omega_m0": omega_m0,
                "omega_r0": omega_r0,
                "error": f"{type(exc).__name__}: {exc}",
                "solver_directory": relative(solver_dir),
            })

    completed = [record for record in assembled if record["status"] == "completed"]
    if not completed:
        raise RuntimeError("no Planck-profile points completed the joint assembly")
    joint_minimum = min(completed, key=lambda record: record["conditional_joint_statistic"])
    planck_minimum = min(completed, key=lambda record: record["planck_minus_2_log_likelihood"])
    late_time_minimum = min(completed, key=lambda record: record["late_time_chi2"]["total"])
    output = {
        "model": "R-Universe covariant KGB",
        "scope": (
            "Conditional joint profile constructed from executed official Planck 2018 likelihood ordinates "
            "and exact-background Pantheon+ plus DESI DR2 BAO plus cosmic-chronometer likelihoods. "
            "All Planck primordial, calibration, and high-ell nuisance inputs are fixed at the values in the source grids; "
            "this is not a posterior, evidence calculation, or complete multi-probe fit."
        ),
        "source_profile_roots": [relative(path) for path in args.profile_roots],
        "source_fixed_inputs": fixed_inputs,
        "late_time_data": {
            "pantheon_plus_count": pantheon["count"],
            "pantheon_plus_table_sha256": pantheon["table_sha256"],
            "pantheon_plus_covariance_sha256": pantheon["covariance_sha256"],
            "desi_dr2_bao_count": bao["count"],
            "desi_dr2_mean_sha256": bao["mean_sha256"],
            "desi_dr2_covariance_sha256": bao["covariance_sha256"],
            "chronometer_count": chronometers["count"],
            "chronometer_sha256": chronometers["sha256"],
            "sn_selection": f"IS_CALIBRATOR=0 and zHD>{args.sn_z_min:g}; intercept analytically marginalized",
        },
        "numerics": {"late_time_distance_integration_nodes": args.integration_nodes},
        "completed_points": len(completed),
        "failed_points": len(assembled) - len(completed),
        "conditional_joint_minimum": joint_minimum,
        "planck_only_minimum_within_assembled_surface": planck_minimum,
        "late_time_only_minimum_within_assembled_surface": late_time_minimum,
        "records": assembled,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "completed_points": output["completed_points"],
        "failed_points": output["failed_points"],
        "conditional_joint_minimum": joint_minimum,
        "planck_only_minimum_within_assembled_surface": planck_minimum,
        "late_time_only_minimum_within_assembled_surface": late_time_minimum,
    }, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
