#!/usr/bin/env python3
"""Compute directly testable RFG-RXi background, tensor, and PPN observables.

This deliberately does not manufacture CMB spectra.  It evaluates only
quantities fixed by the action and the already stated local-GR matching rule:
the homogeneous background, the tensor-amplitude distance ratio, c_T, and
the Cassini likelihood factor.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from multifluid_reduction import Planck2018Reference
from ppn_likelihood import cassini_minus2loglike, ppn_parameters, solar_weyl_ratio
from rfg_regularized import observables
from rfg_xi_completion import XI_AUDIT_VALUES


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "generated" / "tables"


def xi_observable_rows() -> list[dict[str, float]]:
    """Return action-defined observables for the declared Xi benchmarks."""
    reference = Planck2018Reference()
    params = reference.rfg_background_parameters()
    local_ppn = ppn_parameters(solar_weyl_ratio(149597870700.0, reference.h * 100.0))
    scale_factors = np.logspace(-8.0, 0.0, 301)
    rows: list[dict[str, float]] = []
    for xi_completion in XI_AUDIT_VALUES:
        today = observables(1.0, params)
        q_tensor_today = float(today["Q_T"]) + xi_completion
        if q_tensor_today <= 0.0:
            raise RuntimeError("non-positive completed tensor kinetic coefficient today")
        for a in scale_factors:
            background = observables(float(a), params)
            q_tensor = float(background["Q_T"]) + xi_completion
            if q_tensor <= 0.0:
                raise RuntimeError("non-positive completed tensor kinetic coefficient")
            rows.append(
                {
                    "xi_completion": xi_completion,
                    "a": float(a),
                    "z": float(background["z"]),
                    "E": float(background["E"]),
                    "H_km_s_Mpc": reference.h * 100.0 * float(background["E"]),
                    "Omega_m": float(background["Omega_m"]),
                    "Omega_r": float(background["Omega_r"]),
                    "Omega_R": float(background["Omega_R"]),
                    "Q_tensor": q_tensor,
                    "c_T": 1.0,
                    "dL_GW_over_dL_EM": math.sqrt(q_tensor_today / q_tensor),
                    "gamma_PPN": local_ppn["gamma"],
                    "beta_PPN": local_ppn["beta"],
                    "alpha1_PPN": local_ppn["alpha1"],
                    "alpha2_PPN": local_ppn["alpha2"],
                    "Cassini_minus_2loglike": cassini_minus2loglike(local_ppn["gamma"]),
                }
            )
    return rows


def _audit_summary() -> dict[float, dict[str, float]]:
    """Read the exact DAE audit generated earlier in the same validation run."""
    path = TABLES / "rfg_xi_completion_audit.csv"
    if not path.is_file():
        raise RuntimeError("run rfg_xi_completion.py before computing RFG-RXi observables")
    summary: dict[float, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            xi_completion = float(row["xi_completion"])
            current = summary.setdefault(
                xi_completion,
                {
                    "dae_root_bearing_scale_factors": 0.0,
                    "min_DAE_mu_relative": float("inf"),
                    "min_DAE_constraint_determinant": float("inf"),
                },
            )
            current["dae_root_bearing_scale_factors"] = max(
                current["dae_root_bearing_scale_factors"], float(row["root_at_this_a"])
            )
            current["min_DAE_mu_relative"] = min(current["min_DAE_mu_relative"], float(row["mu_zeta_relative"]))
            current["min_DAE_constraint_determinant"] = min(
                current["min_DAE_constraint_determinant"], abs(float(row["constraint_determinant"]))
            )
    return summary


def write_outputs(rows: list[dict[str, float]]) -> tuple[Path, Path]:
    """Write the action-defined observable grid and concise validity summary."""
    TABLES.mkdir(parents=True, exist_ok=True)
    grid_path = TABLES / "rfg_xi_observables.csv"
    with grid_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    audit = _audit_summary()
    summary_path = TABLES / "rfg_xi_validated_factors.csv"
    summary_rows: list[dict[str, float]] = []
    for xi_completion in XI_AUDIT_VALUES:
        branch = [row for row in rows if row["xi_completion"] == xi_completion]
        today = branch[-1]
        z_one = min(branch, key=lambda row: abs(row["z"] - 1.0))
        high_redshift = branch[0]
        if audit[xi_completion]["dae_root_bearing_scale_factors"] != 0.0:
            raise RuntimeError("RFG-RXi DAE audit reported a root-bearing scale factor")
        summary_rows.append(
            {
                "xi_completion": xi_completion,
                "Cassini_minus_2loglike": today["Cassini_minus_2loglike"],
                "c_T_minus_1": today["c_T"] - 1.0,
                "Q_tensor_today": today["Q_tensor"],
                "dL_GW_over_dL_EM_z_approximately_1": z_one["dL_GW_over_dL_EM"],
                "dL_GW_over_dL_EM_high_redshift": high_redshift["dL_GW_over_dL_EM"],
                **audit[xi_completion],
            }
        )
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    return grid_path, summary_path


def main() -> None:
    rows = xi_observable_rows()
    grid_path, summary_path = write_outputs(rows)
    with summary_path.open(newline="", encoding="utf-8") as handle:
        summaries = {
            float(row["xi_completion"]): row for row in csv.DictReader(handle)
        }
    for xi_completion in XI_AUDIT_VALUES:
        summary = summaries[xi_completion]
        print(
            f"Xi={xi_completion:g}: Cassini -2 ln L={float(summary['Cassini_minus_2loglike']):.10f}, "
            f"dL_GW/dL_EM(z~1)={float(summary['dL_GW_over_dL_EM_z_approximately_1']):.10f}"
        )
    print(f"wrote {grid_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
