#!/usr/bin/env python3
"""Profile the Planck absolute-calibration nuisance at fixed spectra.

The cosmological model and its spectra are held fixed.  The script only
evaluates the official Planck likelihood at declared A_planck values, so its
output is a one-dimensional nuisance profile, not a cosmological posterior.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
EVALUATOR = ROOT / "scripts" / "evaluate_planck_2018_fixed.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectra-dir", type=Path, required=True)
    parser.add_argument("--planck-base", type=Path, default=os.environ.get("PLANCK_2018_BASE"))
    parser.add_argument("--clipy-source", type=Path, default=os.environ.get("CLIPY_SOURCE"))
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument(
        "--a-planck-values",
        type=float,
        nargs="+",
        default=[0.994, 0.996, 0.998, 1.000, 1.002, 1.004, 1.006],
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def require_directory(path: Path | None, label: str) -> Path:
    if path is None or not path.is_dir():
        raise FileNotFoundError(f"{label} is not available: {path}")
    return path


def run(command: list[str], cwd: Path, log: Path) -> None:
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    log.write_text(completed.stdout + completed.stderr, encoding="utf-8")
    completed.check_returncode()


def main() -> None:
    args = parse_args()
    if not args.spectra_dir.is_dir():
        raise FileNotFoundError(f"spectrum directory is not available: {args.spectra_dir}")
    args.planck_base = require_directory(args.planck_base, "Planck 2018 distribution")
    args.clipy_source = require_directory(args.clipy_source, "clipy source")
    if not args.python.is_file():
        raise FileNotFoundError(f"Python interpreter is not available: {args.python}")
    if any(value <= 0.0 for value in args.a_planck_values):
        raise ValueError("A_planck values must be positive")

    args.output.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    for value in args.a_planck_values:
        point = args.output / f"a_planck_{value:.6f}.json"
        log = args.output / f"a_planck_{value:.6f}.log"
        command = [
            str(args.python), str(EVALUATOR),
            "--spectra-dir", str(args.spectra_dir),
            "--planck-base", str(args.planck_base),
            "--clipy-source", str(args.clipy_source),
            "--a-planck", format(value, ".17g"),
            "--output", str(point),
        ]
        try:
            run(command, ROOT, log)
            report = json.loads(point.read_text(encoding="utf-8"))
            records.append({
                "status": report["status"],
                "A_planck": value,
                "minus_2_log_likelihood": report["total"]["minus_2_log_likelihood"],
                "components": {
                    name: component["minus_2_log_likelihood"]
                    for name, component in report["components"].items()
                },
            })
        except Exception as exc:
            records.append({
                "status": "failed",
                "A_planck": value,
                "error": f"{type(exc).__name__}: {exc}",
            })

    completed = [record for record in records if record["status"] == "completed"]
    summary = {
        "scope": (
            "One-dimensional A_planck nuisance profile at fixed cosmological spectra; "
            "not a cosmological posterior or a model comparison."
        ),
        "records": records,
        "minimum": min(completed, key=lambda record: record["minus_2_log_likelihood"]) if completed else None,
    }
    summary_path = args.output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
