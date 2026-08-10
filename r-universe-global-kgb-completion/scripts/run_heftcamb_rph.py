#!/usr/bin/env python3
"""Run and convergence-test the native H-EFTCAMB RPH realization of R-Universe.

The source functions are generated from ``ru_kgb.py``.  Each spline resolution
is executed in a separate directory, and the physical CMB and matter outputs
are compared directly.  This is a fixed-parameter solver gate, not a data
likelihood or parameter fit.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "scripts" / "generate_heftcamb_rph.py"
DEFAULT_BINARY = Path("/Users/mpetr/Documents/Codex/2026-08-09/m/work/EFTCAMB/fortran/camb")
DEFAULT_TEMPLATE = Path("/Users/mpetr/Documents/Codex/2026-08-09/m/work/EFTCAMB/fortran/HighLExtrapTemplate_lenspotentialCls.dat")


def read_table(path: Path) -> np.ndarray:
    table = np.loadtxt(path, comments="#")
    if table.ndim != 2 or not np.all(np.isfinite(table)):
        raise RuntimeError(f"invalid numerical output: {path}")
    return table


def fractional_rms(current: np.ndarray, reference: np.ndarray, floor: float) -> tuple[float, float]:
    mask = np.abs(reference) > floor
    if not np.any(mask):
        return 0.0, 0.0
    fractional = (current[mask] - reference[mask]) / reference[mask]
    return float(np.sqrt(np.mean(fractional * fractional))), float(np.max(np.abs(fractional)))


def compare_outputs(reference_dir: Path, candidate_dir: Path) -> dict[str, float]:
    reference_cls = read_table(reference_dir / "ru_kgb_rph_lensedCls.dat")
    candidate_cls = read_table(candidate_dir / "ru_kgb_rph_lensedCls.dat")
    if not np.array_equal(reference_cls[:, 0], candidate_cls[:, 0]):
        raise RuntimeError("CMB multipole grids differ between convergence runs")
    reference_pk = read_table(reference_dir / "ru_kgb_rph_matterpower.dat")
    candidate_pk = read_table(candidate_dir / "ru_kgb_rph_matterpower.dat")
    candidate_pk_values = np.interp(reference_pk[:, 0], candidate_pk[:, 0], candidate_pk[:, 1])
    tt_rms, tt_max = fractional_rms(candidate_cls[:, 1], reference_cls[:, 1], 1.0e-6)
    ee_rms, ee_max = fractional_rms(candidate_cls[:, 2], reference_cls[:, 2], 1.0e-8)
    te_rms, te_max = fractional_rms(candidate_cls[:, 4], reference_cls[:, 4], 1.0e-8)
    pk_rms, pk_max = fractional_rms(candidate_pk_values, reference_pk[:, 1], 1.0e-12)
    return {
        "TT_fractional_rms": tt_rms,
        "TT_fractional_max": tt_max,
        "EE_fractional_rms": ee_rms,
        "EE_fractional_max": ee_max,
        "TE_fractional_rms": te_rms,
        "TE_fractional_max": te_max,
        "Pk_fractional_rms": pk_rms,
        "Pk_fractional_max": pk_max,
    }


def run_one(points: int, binary: Path, template: Path, output_root: Path) -> Path:
    run_dir = output_root / f"nodes_{points}"
    run_dir.mkdir(parents=True, exist_ok=True)
    generation = subprocess.run(
        [sys.executable, str(GENERATOR), "--points", str(points), "--output", str(run_dir)],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    (run_dir / "generation.log").write_text(generation.stdout + generation.stderr, encoding="ascii")
    generation.check_returncode()
    shutil.copyfile(template, run_dir / template.name)
    solver = subprocess.run([str(binary), "ru_kgb_rph.ini"], cwd=run_dir, text=True, capture_output=True)
    (run_dir / "solver.log").write_text(solver.stdout + solver.stderr, encoding="ascii")
    solver.check_returncode()
    required = ["ru_kgb_rph_lensedCls.dat", "ru_kgb_rph_scalCls.dat", "ru_kgb_rph_matterpower.dat"]
    missing = [name for name in required if not (run_dir / name).is_file()]
    if missing:
        raise RuntimeError(f"H-EFTCAMB completed without required outputs: {missing}")
    return run_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=int, nargs="+", default=[201, 401, 601])
    parser.add_argument("--binary", type=Path, default=Path(os.environ.get("HEFTCAMB_BIN", DEFAULT_BINARY)))
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--output", type=Path, default=ROOT / "generated" / "heftcamb" / "convergence")
    args = parser.parse_args()
    nodes = sorted(set(args.nodes))
    if len(nodes) < 2 or any(value < 6 for value in nodes):
        raise ValueError("provide at least two spline resolutions, each >= 6")
    if not args.binary.is_file() or not os.access(args.binary, os.X_OK):
        raise FileNotFoundError(f"H-EFTCAMB executable not available: {args.binary}")
    if not args.template.is_file():
        raise FileNotFoundError(f"H-EFTCAMB template not available: {args.template}")

    runs = {points: run_one(points, args.binary, args.template, args.output) for points in nodes}
    reference_points = nodes[-1]
    comparisons = {
        str(points): compare_outputs(runs[reference_points], runs[points]) for points in nodes[:-1]
    }
    summary = {
        "model": "R-Universe covariant KGB, native H-EFTCAMB RPH execution",
        "scope": "fixed-parameter CMB and linear-matter solver convergence; not a data likelihood or posterior",
        "nodes": nodes,
        "reference_nodes": reference_points,
        "comparisons_to_reference": comparisons,
        "outputs": {str(points): str(directory) for points, directory in runs.items()},
    }
    args.output.mkdir(parents=True, exist_ok=True)
    summary_path = args.output / "heftcamb_convergence.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(json.dumps(summary, indent=2))
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
