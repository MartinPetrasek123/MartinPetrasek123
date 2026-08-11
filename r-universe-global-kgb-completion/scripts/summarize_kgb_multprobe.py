#!/usr/bin/env python3
"""Assemble the executed KGB CMB, late-time, RSD-audit, and action checks.

This script intentionally does not manufacture a global likelihood.  It reads
only executed records and labels the conditional Planck--Pantheon+--DESI--CC
sum separately from the native RSD residual audit and the local PPN gate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"required executed record is missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated-root", type=Path, default=ROOT / "generated")
    parser.add_argument("--output", type=Path, default=ROOT / "generated" / "kgb_multprobe_conditional_summary.json")
    return parser.parse_args()


def minimum(record: dict[str, Any], label: str) -> dict[str, Any]:
    value = record.get("minimum")
    if not isinstance(value, dict) or "minus_2_log_likelihood" not in value:
        raise RuntimeError(f"{label} has no completed likelihood minimum")
    return value


def main() -> None:
    args = parse_args()
    generated = args.generated_root
    selected = generated / "planck_profile_joint_center_refined_0975"

    kgb_calibration = minimum(read_json(selected / "calibration_kgb_fine" / "summary.json"), "KGB calibration")
    lcdm_calibration = minimum(read_json(selected / "calibration_lcdm_fine" / "summary.json"), "LCDM calibration")
    kgb_late = read_json(selected / "kgb_late_time_nodes_32769.json")
    lcdm_late = read_json(selected / "lcdm_late_time.json")
    action = read_json(selected / "best_point_action_validation.json")
    ppn = read_json(selected / "best_point_ppn_screening.json")
    surface = read_json(generated / "joint_planck_late_time_existing_surface.json")
    fine_surface = read_json(generated / "joint_planck_late_time_center_refined_0975.json")
    rsd_primary = read_json(generated / "rsd_native_best_point" / "rsd_native_dz_0p005.json")
    rsd_half_step = read_json(generated / "rsd_native_best_point_dz_0p0025" / "rsd_native_dz_0p0025.json")

    parameters = kgb_late["parameters"]
    lcdm_parameters = lcdm_late["parameters"]
    for key in ("omega_m0", "omega_r0", "omega_R0"):
        if not math.isclose(parameters[key], lcdm_parameters[key], rel_tol=0.0, abs_tol=1.0e-15):
            raise RuntimeError(f"KGB/LCDM late-time inputs disagree in {key}")
    if lcdm_parameters["alpha"] != 0.0:
        raise RuntimeError("matched LCDM record must use alpha=0")

    kgb_planck = float(kgb_calibration["minus_2_log_likelihood"])
    lcdm_planck = float(lcdm_calibration["minus_2_log_likelihood"])
    kgb_late_total = float(kgb_late["chi2"]["total"])
    lcdm_late_total = float(lcdm_late["chi2"]["total"])
    kgb_conditional_total = kgb_planck + kgb_late_total
    lcdm_conditional_total = lcdm_planck + lcdm_late_total

    primary_predictions = rsd_primary["predictions"]
    half_predictions = {row["name"]: row for row in rsd_half_step["predictions"]}
    if set(half_predictions) != {row["name"] for row in primary_predictions}:
        raise RuntimeError("the two RSD finite-difference records do not cover the same observations")
    step_changes = [
        half_predictions[row["name"]]["f_sigma8_native"] - row["f_sigma8_native"] for row in primary_predictions
    ]
    velocity_differences = [
        abs(row["difference"]) for row in rsd_primary["solver"]["f_sigma8_velocity_density_cross_check"]
    ]
    sigma8_differences = [
        abs(row["difference"]) for row in rsd_primary["solver"]["sigma8_power_integral_cross_check"]
    ]

    output = {
        "model": "R-Universe covariant KGB completion",
        "scope": (
            "Executed, fixed-input records only. The Planck--Pantheon+--DESI DR2 BAO--chronometer sum is conditional, "
            "not a posterior, evidence calculation, or complete global fit. Native RSD values are a residual audit and "
            "are intentionally excluded from that sum because the supplied compilation lacks the survey likelihood inputs."
        ),
        "selected_matched_point": {
            **parameters,
            "H0_km_s_Mpc": kgb_late["chi2"]["H0_km_s_Mpc"],
            "r_drag_Mpc": kgb_late["chi2"]["r_drag_Mpc"],
            "rph_spline_nodes": 601,
            "fixed_primordial_inputs": {"omega_b": 0.02237, "A_s": 2.1e-9, "n_s": 0.9649, "tau": 0.0544},
        },
        "action_validation": action,
        "selected_point_ppn_screening_gate": {
            "parameters": ppn["parameters"],
            "r_V_pc": ppn["screening"]["r_V_pc"],
            "absolute_gamma_minus_one_envelope": ppn["PPN_gamma"]["absolute_gamma_minus_one_envelope"],
            "Cassini_absolute_bound": ppn["PPN_gamma"]["Cassini_absolute_bound"],
            "margin_factor": ppn["PPN_gamma"]["margin_factor"],
            "scope": ppn["scope"],
            "source_record": "generated/planck_profile_joint_center_refined_0975/best_point_ppn_screening.json",
        },
        "conditional_planck_late_time_pair": {
            "definition": "Each model uses its own minimum on the executed fixed-spectrum A_planck grid; all remaining cosmological, primordial, recombination, and high-ell nuisance inputs are fixed and matched.",
            "kgb": {
                "A_planck": kgb_calibration["A_planck"],
                "planck_minus_2_log_likelihood": kgb_planck,
                "late_time_chi2": kgb_late["chi2"],
                "conditional_sum": kgb_conditional_total,
            },
            "matched_lcdm": {
                "A_planck": lcdm_calibration["A_planck"],
                "planck_minus_2_log_likelihood": lcdm_planck,
                "late_time_chi2": lcdm_late["chi2"],
                "conditional_sum": lcdm_conditional_total,
            },
            "KGB_minus_LCDM_conditional_sum": kgb_conditional_total - lcdm_conditional_total,
            "source_records": {
                "kgb_calibration": "generated/planck_profile_joint_center_refined_0975/calibration_kgb_fine/summary.json",
                "lcdm_calibration": "generated/planck_profile_joint_center_refined_0975/calibration_lcdm_fine/summary.json",
                "kgb_late_time": "generated/planck_profile_joint_center_refined_0975/kgb_late_time_nodes_32769.json",
                "lcdm_late_time": "generated/planck_profile_joint_center_refined_0975/lcdm_late_time.json",
            },
        },
        "fixed_A_planck_conditional_surfaces": {
            "existing_154_point_surface_minimum": surface["conditional_joint_minimum"],
            "center_refinement_minimum": fine_surface["conditional_joint_minimum"],
            "scope": "A_planck=1 and all other non-displayed inputs fixed; the records establish local behavior only.",
        },
        "native_rsd_audit": {
            "scope": rsd_primary["scope"],
            "record_count": rsd_primary["data"]["record_count"],
            "diag_labelled_count": rsd_primary["data"]["diag_labelled_count"],
            "diagonal_residual_sum_not_a_likelihood": rsd_primary["diagnostic_only"]["diagonal_sum_of_squared_residuals_for_records_labelled_diag"],
            "f_sigma8_range": {
                "minimum": min(row["f_sigma8_native"] for row in primary_predictions),
                "maximum": max(row["f_sigma8_native"] for row in primary_predictions),
            },
            "numerical_checks": {
                "max_abs_sigma8_power_integral_minus_solver_printout": max(sigma8_differences),
                "max_abs_f_sigma8_power_derivative_minus_velocity_density": max(velocity_differences),
                "max_abs_f_sigma8_change_when_dz_is_halved": max(abs(value) for value in step_changes),
                "mean_abs_f_sigma8_change_when_dz_is_halved": sum(abs(value) for value in step_changes) / len(step_changes),
            },
            "source_records": {
                "dz_0p005": "generated/rsd_native_best_point/rsd_native_dz_0p005.json",
                "dz_0p0025": "generated/rsd_native_best_point_dz_0p0025/rsd_native_dz_0p0025.json",
            },
        },
        "empirical_boundary": (
            "No claim of empirical replacement follows from this record. A released full inference must sample the cosmological "
            "and nuisance parameters and include official survey likelihoods, covariance/window/AP treatments, nonlinear modeling, "
            "and a declared comparison protocol."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="ascii")
    print(json.dumps(output, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
