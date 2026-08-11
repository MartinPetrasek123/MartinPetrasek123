#!/usr/bin/env python3
"""Merge conditional KGB joint-profile chunks without recomputing likelihoods."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "generated" / "joint_planck_late_time_existing_surface.json")
    return parser.parse_args()


def point_key(record: dict[str, Any]) -> tuple[float, float, float]:
    return float(record["alpha"]), float(record["omega_m0"]), float(record["omega_r0"])


def main() -> None:
    args = parse_args()
    records_by_point: dict[tuple[float, float, float], dict[str, Any]] = {}
    source_records: list[str] = []
    scopes: list[str] = []
    for path in args.inputs:
        report = json.loads(path.read_text(encoding="utf-8"))
        source_records.append(str(path))
        scopes.append(str(report.get("scope", "")))
        for record in report.get("records", []):
            key = point_key(record)
            old = records_by_point.get(key)
            if old is None or (record.get("status") == "completed" and old.get("status") != "completed"):
                records_by_point[key] = record
    records = [records_by_point[key] for key in sorted(records_by_point)]
    completed = [record for record in records if record.get("status") == "completed"]
    if not completed:
        raise RuntimeError("no completed conditional-joint records")
    output = {
        "model": "R-Universe covariant KGB",
        "scope": scopes[0],
        "source_assemblies": source_records,
        "completed_points": len(completed),
        "failed_points": len(records) - len(completed),
        "conditional_joint_minimum": min(completed, key=lambda row: row["conditional_joint_statistic"]),
        "planck_only_minimum_within_surface": min(completed, key=lambda row: row["planck_minus_2_log_likelihood"]),
        "late_time_only_minimum_within_surface": min(completed, key=lambda row: row["late_time_chi2"]["total"]),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: output[key] for key in (
        "completed_points", "failed_points", "conditional_joint_minimum",
        "planck_only_minimum_within_surface", "late_time_only_minimum_within_surface",
    )}, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
