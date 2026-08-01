#!/usr/bin/env python3
"""Reproduce the MLPerf Tiny v1.3 VIRR_task table from a local results repo.

Usage:
  python reproduce_mlperf_tiny_v1_3_virr.py /path/to/tiny_results_v1.3 out.csv
"""

from __future__ import annotations

import csv
import math
import re
import subprocess
import sys
from pathlib import Path


def commit(repo: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def trial_energies(path: Path) -> list[float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return [float(x) for x in re.findall(r"Energy/Inf\.?\s+:\s+([0-9.]+)\s+uJ/inf\.?", text)]


def median(xs: list[float]) -> float:
    ys = sorted(xs)
    n = len(ys)
    if n == 0:
        return math.nan
    mid = n // 2
    return ys[mid] if n % 2 else 0.5 * (ys[mid - 1] + ys[mid])


def parse_accuracy(path: Path) -> tuple[str, str]:
    if not path.exists():
        return "", ""
    text = path.read_text(encoding="utf-8", errors="replace")
    top = re.search(r"Top 1%\s*=\s*([0-9.]+)", text)
    auc = re.search(r"AUC\s*=\s*([0-9.]+)", text)
    return (top.group(1) if top else "", auc.group(1) if auc else "")


def main() -> None:
    repo = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("work/tiny_results_v1.3")
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("mlperf_tiny_v1_3_virr_reproducibility_table.csv")
    fields = [
        "dataset_version", "commit_hash", "benchmark", "submitter", "system",
        "source_file", "trial_energies_uJ_per_outcome", "median_uJ_per_outcome",
        "VIRR_task_outcomes_per_J", "top1_percent", "auc"
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
        med = median(trials)
        rows.append({
            "dataset_version": "MLPerf Tiny v1.3",
            "commit_hash": commit(repo),
            "benchmark": benchmark,
            "submitter": submitter,
            "system": system,
            "source_file": str(path),
            "trial_energies_uJ_per_outcome": ";".join(f"{x:.8g}" for x in trials),
            "median_uJ_per_outcome": f"{med:.8g}",
            "VIRR_task_outcomes_per_J": f"{1_000_000.0 / med:.8g}",
            "top1_percent": top1,
            "auc": auc,
        })
    rows.sort(key=lambda r: (r["benchmark"], float(r["median_uJ_per_outcome"])))
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fields)
        writer.writeheader()
        writer.writerows(rows)
    print(out)


if __name__ == "__main__":
    main()
