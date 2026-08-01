#!/usr/bin/env python3
"""Reproduce the MLPerf Tiny v1.3 VIRR_task table from a local results repo.

Usage:
  python reproduce_mlperf_tiny_v1_3_virr.py /path/to/tiny_results_v1.3 out.csv
"""

from __future__ import annotations

import csv
import argparse
import math
import re
import subprocess
from pathlib import Path

import numpy as np


DOWNLOAD_DATE = "2026-08-01"

def commit(repo: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def trial_energies(path: Path) -> list[float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return [float(x) for x in re.findall(r"Energy/Inf\.?\s+:\s+([0-9.]+)\s+uJ/inf\.?", text)]


def submitted_median(path: Path, trials: list[float]) -> float:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Median energy cost is\s+([0-9.]+)\s+uJ/inf", text)
    return float(match.group(1)) if match else float(np.median(trials))


def bootstrap_ci_from_trials(trials: list[float], seed: int, n: int = 1200) -> tuple[float, float]:
    if len(trials) < 2:
        return (math.nan, math.nan)
    arr = np.array(trials, dtype=float)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        sample = rng.choice(arr, size=len(arr), replace=True)
        vals.append(1_000_000.0 / np.median(sample))
    return tuple(np.percentile(vals, [2.5, 97.5]))


def parse_accuracy(path: Path) -> tuple[str, str]:
    if not path.exists():
        return "", ""
    text = path.read_text(encoding="utf-8", errors="replace")
    top = re.search(r"Top 1%\s*=\s*([0-9.]+)", text)
    auc = re.search(r"AUC\s*=\s*([0-9.]+)", text)
    return (top.group(1) if top else "", auc.group(1) if auc else "")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reproduce the MLPerf Tiny v1.3 VIRR_task table from a local results checkout."
    )
    parser.add_argument(
        "repo",
        nargs="?",
        default="work/tiny_results_v1.3",
        help="Path to the local MLPerf Tiny v1.3 results repository.",
    )
    parser.add_argument(
        "out",
        nargs="?",
        default="mlperf_tiny_v1_3_virr_reproducibility_table.csv",
        help="Output CSV path.",
    )
    parser.add_argument(
        "--source-prefix",
        default=None,
        help=(
            "Prefix used in the source_file column. By default the input repo path is used. "
            "For byte-identical S1 reproduction from this supplement, run from the manuscript "
            "working directory with repo=work/tiny_results_v1.3."
        ),
    )
    args = parser.parse_args()
    repo = Path(args.repo)
    out = Path(args.out)
    fields = [
        "dataset_version", "download_date", "commit_hash", "benchmark", "submitter", "system",
        "source_file", "trial_energies_uJ_per_inf", "median_uJ_per_inf",
        "VIRR_task_inf_per_J", "bootstrap_95pct_inf_per_J", "quality_statistic"
    ]
    rows = []
    for path in sorted(repo.glob("closed/*/results/*/*/energy/results.txt")):
        trials = trial_energies(path)
        if not trials:
            continue
        parts = path.parts
        idx = parts.index("closed")
        submitter, system, benchmark = parts[idx + 1], parts[idx + 3], parts[idx + 4]
        top1, auc = parse_accuracy(path.parents[1] / "accuracy" / "results.txt")
        med = submitted_median(path, trials)
        virr = 1_000_000.0 / med
        source_prefix = Path(args.source_prefix) if args.source_prefix else repo
        source_file = str(source_prefix / path.relative_to(repo))
        q = ""
        if top1:
            q = f"Top-1={top1}"
        if auc:
            q = (q + "; " if q else "") + f"AUC={auc}"
        rows.append({
            "dataset_version": "MLPerf Tiny v1.3",
            "download_date": DOWNLOAD_DATE,
            "commit_hash": commit(repo),
            "benchmark": benchmark,
            "submitter": submitter,
            "system": system,
            "source_file": source_file,
            "trial_energies_uJ_per_inf": ";".join(f"{x:.8g}" for x in trials),
            "median_uJ_per_inf": f"{med:.8g}",
            "VIRR_task_inf_per_J": f"{virr:.8g}",
            "quality_statistic": q,
            "_virr": virr,
        })
    rows.sort(key=lambda r: r["_virr"], reverse=True)
    for i, row in enumerate(rows, 1):
        trials = [float(x) for x in row["trial_energies_uJ_per_inf"].split(";") if x]
        lo, hi = bootstrap_ci_from_trials(trials, seed=200 + i, n=1200)
        row["bootstrap_95pct_inf_per_J"] = (
            "not estimable from a single submitted trial"
            if len(trials) < 2
            else f"{lo:.8g}-{hi:.8g}"
        )
        del row["_virr"]
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fields)
        writer.writeheader()
        writer.writerows(rows)
    print(out)


if __name__ == "__main__":
    main()
