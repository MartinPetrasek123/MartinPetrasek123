#!/usr/bin/env python3
"""Full multiplier-sector kinetic gate for locked R-Universe.

The integrated-by-parts multiplier action contains

    - Z(R)/2 (grad R)^2 - ell^2 grad(lambda).grad(R)

For homogeneous perturbations around any background, the kinetic Hessian in
the field velocities (Rdot, lambdadot) is

    K = [[Z(R), ell^2],
         [ell^2, 0]]

up to the common positive volume factor a^3.  Its determinant is -ell^4, so
one eigenvalue is negative for every finite nonzero ell.  This script evaluates
the gate on the previously integrated locked background and records the
background-independent conclusion.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BACKGROUND = ROOT / "tables" / "locked_r_universe_background.csv"
RESULTS = ROOT / "code" / "full_lambda_stability_gate_results.json"
CSV = ROOT / "tables" / "full_lambda_stability_gate.csv"
FIG_PDF = ROOT / "figures" / "fig02_full_lambda_ghost_gate.pdf"
FIG_PNG = ROOT / "figures" / "fig02_full_lambda_ghost_gate.png"


def kinetic_eigenvalues(z: np.ndarray, ell: float) -> tuple[np.ndarray, np.ndarray]:
    root = np.sqrt(z * z + 4.0 * ell**4)
    return 0.5 * (z + root), 0.5 * (z - root)


def main() -> None:
    df = pd.read_csv(BACKGROUND)
    # The current locked branch selected ell=0.7; the sign result is independent
    # of this value as long as ell is finite and nonzero.
    ell = 0.7
    z = df["Z"].to_numpy()
    kp, km = kinetic_eigenvalues(z, ell)
    out_df = pd.DataFrame(
        {
            "a": df["a"],
            "Z": z,
            "ell": ell,
            "K_plus": kp,
            "K_minus": km,
            "det_K": -ell**4,
        }
    )
    out_df.to_csv(CSV, index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(out_df["a"], out_df["K_plus"], label=r"$K_+$")
    ax.plot(out_df["a"], out_df["K_minus"], label=r"$K_-$")
    ax.axhline(0.0, color="0.2", lw=0.9)
    ax.set_xscale("log")
    ax.set_xlabel("a")
    ax.set_ylabel("kinetic eigenvalue")
    ax.set_title("full multiplier kinetic gate")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_PDF)
    fig.savefig(FIG_PNG, dpi=180)
    plt.close(fig)

    result = {
        "gate": "full multiplier scalar kinetic Hessian",
        "kinetic_matrix": "[[Z(R), ell^2], [ell^2, 0]]",
        "determinant": "-ell^4",
        "background_independent_conclusion": (
            "For any finite nonzero ell and positive Z, det(K)<0, so the two "
            "eigenvalues have opposite signs. The Lagrange-multiplier action "
            "contains a ghost scalar degree of freedom before data fitting."
        ),
        "current_branch": {
            "ell": ell,
            "min_Z": float(np.min(z)),
            "min_K_plus": float(np.min(kp)),
            "max_K_minus": float(np.max(km)),
            "min_K_minus": float(np.min(km)),
            "all_K_plus_positive": bool(np.all(kp > 0)),
            "all_K_minus_negative": bool(np.all(km < 0)),
        },
        "verdict": "FAIL: locked Lagrange-multiplier finite-window action is not viable in this form",
        "possible_next_theory_move": (
            "Replace the hard multiplier by a healthy auxiliary field or a "
            "degenerate constrained construction whose kinetic matrix has no "
            "negative propagating eigenvalue."
        ),
        "outputs": {
            "csv": str(CSV.relative_to(ROOT)),
            "figure": str(FIG_PDF.relative_to(ROOT)),
        },
    }
    RESULTS.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
