#!/usr/bin/env python3
"""Run a declared KGB Planck profile grid and a matched fixed LCDM reference.

This script is deliberately not an MCMC or an evidence calculation.  It holds
the primordial sector and Planck absolute calibration fixed, evaluates a small
two-dimensional grid in the R-Universe action parameters (alpha, Omega_m0),
and reports the official Planck 2018 likelihood at every completed point.  A
standard LCDM spectrum is evaluated at every sampled Omega_m0 with the same
physical matter, radiation, primordial, and calibration inputs.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

from ru_kgb import RUKGBParams


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "scripts" / "generate_heftcamb_rph.py"
EVALUATOR = ROOT / "scripts" / "evaluate_planck_2018_fixed.py"
DEFAULT_BINARY = Path("/Users/mpetr/Documents/Codex/2026-08-09/m/work/EFTCAMB/fortran/camb")
DEFAULT_TEMPLATE = Path(
    "/Users/mpetr/Documents/Codex/2026-08-09/m/work/EFTCAMB/fortran/HighLExtrapTemplate_lenspotentialCls.dat"
)


def parse_args() -> argparse.Namespace:
    defaults = RUKGBParams()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path(os.environ.get("HEFTCAMB_BIN", DEFAULT_BINARY)))
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--planck-base", type=Path, default=os.environ.get("PLANCK_2018_BASE"))
    parser.add_argument("--clipy-source", type=Path, default=os.environ.get("CLIPY_SOURCE"))
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--a-planck", type=float, default=1.0)
    parser.add_argument("--points", type=int, default=601)
    parser.add_argument(
        "--alpha-values",
        type=float,
        nargs="+",
        default=[0.45, defaults.alpha, 0.55],
    )
    parser.add_argument(
        "--omega-m0-values",
        type=float,
        nargs="+",
        default=[0.28, defaults.omega_m0, 0.325],
    )
    parser.add_argument("--omega-r0", type=float, default=defaults.omega_r0)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse completed planck.json files in --output instead of rerunning them.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "generated" / "planck_profile_grid",
    )
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not available: {path}")


def require_directory(path: Path | None, label: str) -> Path:
    if path is None or not path.is_dir():
        raise FileNotFoundError(f"{label} is not available: {path}")
    return path


def run(command: list[str], cwd: Path, log: Path) -> None:
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    log.write_text(completed.stdout + completed.stderr, encoding="utf-8")
    completed.check_returncode()


def solver_outputs(directory: Path) -> bool:
    return all((directory / name).is_file() for name in (
        "ru_kgb_rph_lensedCls.dat",
        "ru_kgb_rph_lenspotentialCls.dat",
        "ru_kgb_rph_scalCls.dat",
    ))


def make_kgb_input(run_dir: Path, args: argparse.Namespace, params: RUKGBParams) -> None:
    command = [
        str(args.python), str(GENERATOR),
        "--points", str(args.points),
        "--omega-m0", format(params.omega_m0, ".17g"),
        "--omega-r0", format(params.omega_r0, ".17g"),
        "--alpha", format(params.alpha, ".17g"),
        "--output", str(run_dir),
    ]
    run(command, ROOT, run_dir / "generation.log")
    shutil.copyfile(args.template, run_dir / args.template.name)


def make_lcdm_input(run_dir: Path, args: argparse.Namespace, params: RUKGBParams) -> None:
    """Use the common physical input block, removing all RPH/EFT directives."""
    make_kgb_input(run_dir, args, params)
    source = run_dir / "ru_kgb_rph.ini"
    remove_prefixes = ("EFT", "AltParEFT", "RPH", "model_background")
    retained: list[str] = []
    for line in source.read_text(encoding="ascii").splitlines():
        key = line.split("=", 1)[0].strip()
        if key.startswith(remove_prefixes):
            continue
        retained.append(line)
    retained.extend((
        "# Matched standard reference: no modified-gravity EFT operators.",
        "EFTflag = 0",
    ))
    source.write_text("\n".join(retained) + "\n", encoding="ascii")


def evaluate_point(run_dir: Path, args: argparse.Namespace) -> dict[str, Any]:
    run([str(args.binary), "ru_kgb_rph.ini"], run_dir, run_dir / "solver.log")
    if not solver_outputs(run_dir):
        raise RuntimeError(f"solver did not emit the required spectra in {run_dir}")
    result = run_dir / "planck.json"
    command = [
        str(args.python), str(EVALUATOR),
        "--spectra-dir", str(run_dir),
        "--planck-base", str(args.planck_base),
        "--clipy-source", str(args.clipy_source),
        "--a-planck", format(args.a_planck, ".17g"),
        "--output", str(result),
    ]
    run(command, ROOT, run_dir / "planck.log")
    return json.loads(result.read_text(encoding="utf-8"))


def point_tag(params: RUKGBParams) -> str:
    return f"kgb_alpha_{params.alpha:.8f}_omega_m0_{params.omega_m0:.8f}"


def completed_record(params: RUKGBParams, report: dict[str, Any]) -> dict[str, Any]:
    return {
        "status": report["status"],
        "alpha": params.alpha,
        "omega_m0": params.omega_m0,
        "omega_r0": params.omega_r0,
        "omega_R0": params.omega_R0,
        "minus_2_log_likelihood": report["total"]["minus_2_log_likelihood"],
        "components": {
            key: value["minus_2_log_likelihood"]
            for key, value in report["components"].items()
        },
    }


def main() -> None:
    args = parse_args()
    require_file(args.binary, "H-EFTCAMB binary")
    require_file(args.template, "H-EFTCAMB high-ell template")
    require_file(args.python, "Python interpreter")
    args.planck_base = require_directory(args.planck_base, "Planck 2018 distribution")
    args.clipy_source = require_directory(args.clipy_source, "clipy source")
    if args.points < 6:
        raise ValueError("at least six RPH spline nodes are required")

    args.output.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    for alpha in args.alpha_values:
        for omega_m0 in args.omega_m0_values:
            params = RUKGBParams(omega_m0=omega_m0, omega_r0=args.omega_r0, alpha=alpha)
            params.validate()
            run_dir = args.output / point_tag(params)
            run_dir.mkdir(parents=True, exist_ok=True)
            try:
                cached = run_dir / "planck.json"
                if args.resume and cached.is_file():
                    report = json.loads(cached.read_text(encoding="utf-8"))
                    if report.get("status") != "completed":
                        raise RuntimeError(f"cached likelihood did not complete: {cached}")
                else:
                    make_kgb_input(run_dir, args, params)
                    report = evaluate_point(run_dir, args)
                records.append(completed_record(params, report))
            except Exception as exc:  # Record a rejected point rather than silently dropping it.
                records.append({
                    "status": "failed",
                    "alpha": params.alpha,
                    "omega_m0": params.omega_m0,
                    "omega_r0": params.omega_r0,
                    "omega_R0": params.omega_R0,
                    "error": f"{type(exc).__name__}: {exc}",
                })

    lcdm_by_omega: dict[str, dict[str, Any]] = {}
    for omega_m0 in sorted(set(args.omega_m0_values)):
        params = RUKGBParams(omega_m0=omega_m0, omega_r0=args.omega_r0, alpha=0.0)
        key = format(params.omega_m0, ".17g")
        lcdm_dir = args.output / f"lcdm_omega_m0_{params.omega_m0:.8f}"
        lcdm_dir.mkdir(parents=True, exist_ok=True)
        try:
            cached = lcdm_dir / "planck.json"
            if args.resume and cached.is_file():
                lcdm_report = json.loads(cached.read_text(encoding="utf-8"))
                if lcdm_report.get("status") != "completed":
                    raise RuntimeError(f"cached likelihood did not complete: {cached}")
            else:
                make_lcdm_input(lcdm_dir, args, params)
                lcdm_report = evaluate_point(lcdm_dir, args)
            lcdm_by_omega[key] = {
                "status": lcdm_report["status"],
                "omega_m0": params.omega_m0,
                "omega_r0": params.omega_r0,
                "minus_2_log_likelihood": lcdm_report["total"]["minus_2_log_likelihood"],
                "components": {
                    name: value["minus_2_log_likelihood"]
                    for name, value in lcdm_report["components"].items()
                },
            }
        except Exception as exc:
            lcdm_by_omega[key] = {
                "status": "failed",
                "omega_m0": params.omega_m0,
                "omega_r0": params.omega_r0,
                "error": f"{type(exc).__name__}: {exc}",
            }

    for record in records:
        reference = lcdm_by_omega.get(format(record["omega_m0"], ".17g"))
        if record["status"] == "completed" and reference is not None and reference["status"] == "completed":
            record["delta_minus_2_log_likelihood_vs_matched_lcdm"] = (
                record["minus_2_log_likelihood"] - reference["minus_2_log_likelihood"]
            )

    completed = [record for record in records if record["status"] == "completed"]
    minimum = min(completed, key=lambda record: record["minus_2_log_likelihood"]) if completed else None
    summary: dict[str, Any] = {
        "model": "R-Universe covariant KGB",
        "scope": (
            "Declared two-parameter Planck profile grid with fixed primordial and calibration inputs; "
            "not a posterior, Bayesian evidence calculation, or full model comparison."
        ),
        "fixed_inputs": {
            "omega_b": 0.02237,
            "A_s": 2.1e-9,
            "n_s": 0.9649,
            "tau": 0.0544,
            "A_planck": args.a_planck,
            "rph_spline_nodes": args.points,
        },
        "grid": {
            "alpha_values": args.alpha_values,
            "omega_m0_values": args.omega_m0_values,
            "omega_r0": args.omega_r0,
            "records": records,
            "completed_points": len(completed),
            "grid_minimum": minimum,
        },
        "lcdm_matched_fixed_references": list(lcdm_by_omega.values()),
    }
    if minimum is not None:
        reference = lcdm_by_omega.get(format(minimum["omega_m0"], ".17g"))
        if reference is not None and reference["status"] == "completed":
            summary["grid_minimum_minus_matched_lcdm_fixed_delta_minus_2_log_likelihood"] = (
                minimum["minus_2_log_likelihood"] - reference["minus_2_log_likelihood"]
            )
    result = args.output / "summary.json"
    result.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"wrote {result}")


if __name__ == "__main__":
    main()
