#!/usr/bin/env python3
"""Run an extended action-level DAE root audit for the RFG-RXi benchmarks."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from multifluid_reduction import Planck2018Reference
from rfg_dae_closure import curvature_constraint, quadratic_blocks_with_time_derivative
from rfg_xi_completion import XI_AUDIT_VALUES


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "generated" / "tables"


def _relative_mu(closure: dict[str, object]) -> float:
    mass = np.asarray(closure["M"], dtype=float)
    mixing_dot = np.asarray(closure["B_dot"], dtype=float)
    return abs(float(closure["mu_zeta"])) / max(1.0, abs(mass[0, 0]), abs(mixing_dot[0, 0]))


def _root_on_k_grid(
    a: float,
    k_values: np.ndarray,
    *,
    params: object,
    species: tuple[object, ...],
    xi_completion: float,
) -> tuple[float | None, float]:
    """Find every sign-changing root bracket and retain the closest approach."""
    left_k = float(k_values[0])
    left = curvature_constraint(a, left_k, params, species, xi_completion=xi_completion)
    left_mu = float(left["mu_zeta"])
    min_relative = _relative_mu(left)
    for current_k_value in k_values[1:]:
        right_k = float(current_k_value)
        right = curvature_constraint(a, right_k, params, species, xi_completion=xi_completion)
        right_mu = float(right["mu_zeta"])
        min_relative = min(min_relative, _relative_mu(right))
        if left_mu == 0.0:
            return left_k, min_relative
        if left_mu * right_mu < 0.0:
            for _ in range(96):
                middle_k = math.sqrt(left_k * right_k)
                middle = curvature_constraint(a, middle_k, params, species, xi_completion=xi_completion)
                middle_mu = float(middle["mu_zeta"])
                if left_mu * middle_mu <= 0.0:
                    right_k, right_mu = middle_k, middle_mu
                else:
                    left_k, left_mu = middle_k, middle_mu
            return math.sqrt(left_k * right_k), min_relative
        left_k, left_mu = right_k, right_mu
    return None, min_relative


def main() -> None:
    reference = Planck2018Reference()
    params = reference.rfg_background_parameters()
    species = reference.species()
    a_values = np.logspace(-8.0, 0.0, 193)
    k_values = np.logspace(-5.0, 7.0, 601)
    rows: list[dict[str, float | str]] = []
    for xi_completion in XI_AUDIT_VALUES:
        for a in a_values:
            blocks = quadratic_blocks_with_time_derivative(
                float(a), 1.0, params, species, xi_completion=xi_completion
            )
            root, min_relative = _root_on_k_grid(
                float(a),
                k_values,
                params=params,
                species=species,
                xi_completion=xi_completion,
            )
            min_q_tensor = float(blocks["rate"]["Q"]) + xi_completion
            min_determinant = abs(float(np.linalg.det(np.asarray(blocks["A"], dtype=float)[:2, :2])))
            rows.append(
                {
                    "xi_completion": xi_completion,
                    "a": float(a),
                    "z": 1.0 / float(a) - 1.0,
                    "root_k_over_H0": "" if root is None else root,
                    "min_mu_zeta_relative_on_k_grid": min_relative,
                    "min_Q_tensor_on_k_grid": min_q_tensor,
                    "min_constraint_determinant_on_k_grid": min_determinant,
                }
            )
    output = TABLES / "rfg_xi_dense_dae_audit.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    roots = [row for row in rows if row["root_k_over_H0"] != ""]
    if roots:
        raise RuntimeError(f"RFG-RXi dense audit found {len(roots)} curvature-constraint roots")
    if min(float(row["min_Q_tensor_on_k_grid"]) for row in rows) <= 0.0:
        raise RuntimeError("RFG-RXi dense audit found non-positive tensor kinetic coefficient")
    if min(float(row["min_constraint_determinant_on_k_grid"]) for row in rows) <= 0.0:
        raise RuntimeError("RFG-RXi dense audit found singular lapse-shift block")
    print(f"RFG-RXi dense DAE points               = {len(rows)} x {len(k_values)}")
    print(f"sampled curvature-constraint roots      = {len(roots)}")
    print(f"min normalized |mu_zeta|                = {min(float(row['min_mu_zeta_relative_on_k_grid']) for row in rows):.12e}")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
