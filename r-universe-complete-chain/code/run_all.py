#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable

STEPS = [
    "extended_fit.py",
    "des_dovekie_fit.py",
    "derived_predictions.py",
    "numerical_validation.py",
    "boltzmann_camb.py",
    "robustness_suite.py",
    "profile_likelihood.py",
    "make_figures_and_tables.py",
]


def main():
    start = time.time()
    for step in STEPS:
        t0 = time.time()
        print(f"[run_all] {step}")
        subprocess.run([PY, str(ROOT / "code" / step)], cwd=ROOT, check=True)
        print(f"[run_all] {step} completed in {time.time() - t0:.2f} s")
    print(f"[run_all] all steps completed in {time.time() - start:.2f} s")


if __name__ == "__main__":
    main()
