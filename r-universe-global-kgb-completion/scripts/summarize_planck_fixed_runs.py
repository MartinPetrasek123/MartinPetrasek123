#!/usr/bin/env python3
"""Collect fixed-point Planck likelihood runs and their numerical checks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "generated"


def read_report(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"Required Planck run is missing: {path}")
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("status") != "completed":
        raise RuntimeError(f"Planck run did not complete: {path}")
    return report


def chi2(report: dict) -> float:
    return float(report["total"]["minus_2_log_likelihood"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=GENERATED / "planck_2018_fixed_summary.json"
    )
    args = parser.parse_args()

    reports = {
        "nodes_201": read_report(GENERATED / "planck_2018_fixed_loglike_nodes_201.json"),
        "nodes_601": read_report(GENERATED / "planck_2018_fixed_loglike_nodes_601.json"),
        "turn_on_a_1e-4": read_report(GENERATED / "planck_2018_turn_on_a_1e-4.json"),
        "turn_on_a_1e-3": read_report(GENERATED / "planck_2018_turn_on_a_1e-3.json"),
    }
    reference = reports["nodes_601"]
    summary = {
        "status": "completed",
        "scope": reference["scope"],
        "reference_fixed_point": {
            "spline_nodes": 601,
            "a_planck_fixed": reference["inputs"]["a_planck_fixed"],
            "components": reference["components"],
            "total": reference["total"],
        },
        "numerical_checks": {
            "spline_201_minus_601_delta_minus_2_log_likelihood": (
                chi2(reports["nodes_201"]) - chi2(reference)
            ),
            "turn_on_1e-3_minus_1e-4_delta_minus_2_log_likelihood": (
                chi2(reports["turn_on_a_1e-3"]) - chi2(reports["turn_on_a_1e-4"])
            ),
            "turn_on_values": [1.0e-4, 1.0e-3],
        },
        "provenance": {
            "evaluator": "scripts/evaluate_planck_2018_fixed.py",
            "solver_generator": "scripts/generate_heftcamb_rph.py",
            "solver_convergence": "generated/heftcamb/convergence/heftcamb_convergence.json",
            "reports": {
                name: report["inputs"]["spectra_dir"] for name, report in reports.items()
            },
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
