#!/usr/bin/env python3
"""Evaluate one sampled covariant-KGB Planck plus late-time likelihood point.

This is the expensive, action-faithful unit used by a posterior sampler.  It
creates the exact RPH functions for the supplied cosmological and primordial
parameters, executes H-EFTCAMB, evaluates the official Planck 2018 likelihood
objects, and evaluates the exact-background Pantheon+, DESI DR2 BAO, and
chronometer terms.  It deliberately excludes the incomplete RSD compilation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "scripts" / "generate_heftcamb_rph.py"
PLANCK_EVALUATOR = ROOT / "scripts" / "evaluate_planck_2018_fixed.py"
LATE_TIME_EVALUATOR = ROOT / "scripts" / "evaluate_kgb_late_time.py"
DEFAULT_BINARY = Path("/Users/mpetr/Documents/Codex/2026-08-09/m/work/EFTCAMB/fortran/camb")
DEFAULT_TEMPLATE = Path("/Users/mpetr/Documents/Codex/2026-08-09/m/work/EFTCAMB/fortran/HighLExtrapTemplate_lenspotentialCls.dat")


def finite_positive(value: float, name: str) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="Immutable directory for this one parameter point.")
    parser.add_argument("--binary", type=Path, default=Path(os.environ.get("HEFTCAMB_BIN", DEFAULT_BINARY)))
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--clipy-source", type=Path, required=True)
    parser.add_argument("--planck-base", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(sys.executable), help="Python interpreter with clipy and scipy.")
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--omega-m0", type=float, required=True)
    parser.add_argument("--omega-r0", type=float, required=True)
    parser.add_argument("--ombh2", type=float, required=True)
    parser.add_argument("--scalar-amp", type=float, required=True)
    parser.add_argument("--scalar-spectral-index", type=float, required=True)
    parser.add_argument("--tau", type=float, required=True)
    parser.add_argument("--a-planck", type=float, required=True)
    parser.add_argument("--spline-points", type=int, default=601)
    parser.add_argument("--late-time-integration-nodes", type=int, default=16385)
    parser.add_argument("--reuse-completed", action="store_true", help="Return a matching completed summary without rerunning.")
    parser.add_argument(
        "--compact-artifacts",
        action="store_true",
        help="Losslessly archive solver artifacts after a completed likelihood point.",
    )
    return parser.parse_args()


def run_logged(command: list[str], cwd: Path, log_path: Path) -> None:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    result.check_returncode()


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise RuntimeError(f"expected an object in {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def declared_inputs(args: argparse.Namespace) -> dict[str, float | int]:
    return {
        "alpha": args.alpha,
        "omega_m0": args.omega_m0,
        "omega_r0": args.omega_r0,
        "ombh2": args.ombh2,
        "scalar_amp": args.scalar_amp,
        "scalar_spectral_index": args.scalar_spectral_index,
        "tau": args.tau,
        "a_planck": args.a_planck,
        "spline_points": args.spline_points,
        "late_time_integration_nodes": args.late_time_integration_nodes,
    }


def validate_args(args: argparse.Namespace) -> None:
    for path, label in (
        (args.binary, "H-EFTCAMB binary"),
        (args.template, "high-ell extrapolation template"),
        (args.python, "Python interpreter"),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{label} is unavailable: {path}")
    if not os.access(args.binary, os.X_OK):
        raise PermissionError(f"H-EFTCAMB binary is not executable: {args.binary}")
    for path, label in (
        (args.clipy_source, "clipy source"),
        (args.planck_base, "Planck 2018 likelihood data"),
        (args.data_root, "late-time data root"),
    ):
        if not path.is_dir():
            raise FileNotFoundError(f"{label} is unavailable: {path}")
    if not (0.0 <= args.alpha < 1.0):
        raise ValueError("alpha must lie in [0, 1)")
    finite_positive(args.omega_m0, "omega-m0")
    finite_positive(args.omega_r0, "omega-r0")
    if args.omega_m0 + args.omega_r0 >= 1.0:
        raise ValueError("omega-m0 plus omega-r0 must be below one")
    finite_positive(args.ombh2, "ombh2")
    finite_positive(args.scalar_amp, "scalar-amp")
    if not (0.0 < args.scalar_spectral_index < 2.0):
        raise ValueError("scalar-spectral-index must lie in (0, 2)")
    if not (0.0 < args.tau < 1.0):
        raise ValueError("tau must lie in (0, 1)")
    finite_positive(args.a_planck, "a-planck")
    if args.spline_points < 201:
        raise ValueError("spline-points must be at least 201 for a likelihood evaluation")
    if args.late_time_integration_nodes < 1025:
        raise ValueError("late-time-integration-nodes must be at least 1025")


def main() -> None:
    args = parse_args()
    validate_args(args)
    inputs = declared_inputs(args)
    output_dir = args.output_dir.resolve()
    summary_path = output_dir / "joint_point.json"
    if summary_path.is_file():
        existing = read_json(summary_path)
        if existing.get("inputs") != inputs:
            raise RuntimeError(f"output directory already belongs to another parameter point: {output_dir}")
        if args.reuse_completed and existing.get("status") == "completed":
            print(json.dumps(existing, indent=2))
            return
        raise RuntimeError(f"completed point already exists; pass --reuse-completed: {output_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    ini_path = output_dir / "ru_kgb_rph.ini"
    try:
        run_logged(
            [
                str(args.python),
                str(GENERATOR),
                "--points", str(args.spline_points),
                "--omega-m0", format(args.omega_m0, ".17g"),
                "--omega-r0", format(args.omega_r0, ".17g"),
                "--alpha", format(args.alpha, ".17g"),
                "--ombh2", format(args.ombh2, ".17g"),
                "--scalar-amp", format(args.scalar_amp, ".17g"),
                "--scalar-spectral-index", format(args.scalar_spectral_index, ".17g"),
                "--tau", format(args.tau, ".17g"),
                "--output", str(output_dir),
            ],
            ROOT,
            output_dir / "generation.log",
        )
        if not ini_path.is_file():
            raise RuntimeError("generator completed without an H-EFTCAMB input file")
        shutil.copyfile(args.template, output_dir / args.template.name)
        run_logged([str(args.binary), ini_path.name], output_dir, output_dir / "solver.log")
        required_outputs = ("ru_kgb_rph_lensedCls.dat", "ru_kgb_rph_lenspotentialCls.dat", "ru_kgb_rph_matterpower.dat")
        missing = [name for name in required_outputs if not (output_dir / name).is_file()]
        if missing:
            raise RuntimeError(f"H-EFTCAMB completed without required outputs: {missing}")
        planck_path = output_dir / "planck.json"
        run_logged(
            [
                str(args.python), str(PLANCK_EVALUATOR),
                "--spectra-dir", str(output_dir),
                "--planck-base", str(args.planck_base),
                "--clipy-source", str(args.clipy_source),
                "--a-planck", format(args.a_planck, ".17g"),
                "--output", str(planck_path),
            ],
            ROOT,
            output_dir / "planck.log",
        )
        late_path = output_dir / "late_time.json"
        run_logged(
            [
                str(args.python), str(LATE_TIME_EVALUATOR),
                "--data-root", str(args.data_root),
                "--solver-dir", str(output_dir),
                "--omega-m0", format(args.omega_m0, ".17g"),
                "--omega-r0", format(args.omega_r0, ".17g"),
                "--alpha", format(args.alpha, ".17g"),
                "--integration-nodes", str(args.late_time_integration_nodes),
                "--output", str(late_path),
            ],
            ROOT,
            output_dir / "late_time.log",
        )
        planck = read_json(planck_path)
        late_time = read_json(late_path)
        planck_value = float(planck["total"]["minus_2_log_likelihood"])
        late_value = float(late_time["chi2"]["total"])
        if not (math.isfinite(planck_value) and math.isfinite(late_value)):
            raise RuntimeError("one or more likelihood terms are non-finite")
        summary: dict[str, Any] = {
            "status": "completed",
            "model": "R-Universe covariant KGB",
            "scope": (
                "Exact one-point Planck 2018 Plik-lite plus low-ell plus lensing, including the documented "
                "A_planck=1.0000+-0.0025 calibration prior, and Pantheon+ plus DESI DR2 BAO plus chronometer "
                "likelihood ordinate. This is the action-faithful evaluator for a posterior sampler, not itself a "
                "posterior, evidence, RSD likelihood, or full survey analysis."
            ),
            "inputs": inputs,
            "paths": {
                "ini": ini_path.name,
                "solver_log": "solver.log",
                "planck": planck_path.name,
                "late_time": late_path.name,
            },
            "execution_provenance": {
                "generator_sha256": sha256(GENERATOR),
                "planck_wrapper_sha256": sha256(PLANCK_EVALUATOR),
                "late_time_evaluator_sha256": sha256(LATE_TIME_EVALUATOR),
                "heftcamb_binary_sha256": sha256(args.binary),
                "heftcamb_template_sha256": sha256(args.template),
            },
            "minus_2_log_likelihood": {
                "planck": planck_value,
                "late_time": late_value,
                "total": planck_value + late_value,
            },
        }
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        if args.compact_artifacts:
            try:
                if str(ROOT) not in sys.path:
                    sys.path.insert(0, str(ROOT))
                from inference.archive_kgb_joint_cache import compact

                compact(
                    output_dir,
                    prune=True,
                    profile=os.environ.get("KGB_ARCHIVE_PROFILE", "full"),
                )
            except Exception as archive_error:
                # Storage compaction is never allowed to alter a valid datum.
                (output_dir / "archive_failure.json").write_text(
                    json.dumps({"error": f"{type(archive_error).__name__}: {archive_error}"}, indent=2) + "\n",
                    encoding="utf-8",
                )
        print(json.dumps(summary, indent=2))
    except Exception as error:
        failure = {
            "status": "failed",
            "model": "R-Universe covariant KGB",
            "inputs": inputs,
            "error": f"{type(error).__name__}: {error}",
        }
        (output_dir / "joint_point_failure.json").write_text(json.dumps(failure, indent=2) + "\n", encoding="utf-8")
        raise


if __name__ == "__main__":
    main()
