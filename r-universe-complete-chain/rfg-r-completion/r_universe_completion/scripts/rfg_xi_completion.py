#!/usr/bin/env python3
"""Audit the background-preserving Xi(R3+sigma^2) R-Universe completion."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from multifluid_reduction import Planck2018Reference, normalized_inertia
from rfg_dae_closure import curvature_constraint, find_mu_root, quadratic_blocks_with_time_derivative


ROOT = Path(__file__).resolve().parents[1]
XI_REFERENCE = 1.0
XI_CROSSCHECK = 2.0
XI_AUDIT_VALUES = (XI_REFERENCE, XI_CROSSCHECK)


def main() -> None:
    reference = Planck2018Reference()
    params = reference.rfg_background_parameters()
    species = reference.species()
    rows: list[dict[str, float | int | str]] = []
    maximum_kinetic_null = 0.0
    maximum_antisymmetry = 0.0
    minimum_q_tensor = float("inf")
    minimum_constraint_determinant = float("inf")
    minimum_mu_relative = float("inf")

    root_counts: dict[float, int] = {}
    for xi_completion in XI_AUDIT_VALUES:
        root_count = 0
        for a in np.logspace(-7.0, 0.0, 49):
            root = find_mu_root(float(a), params, species, xi_completion=xi_completion)
            root_count += int(root is not None)
            for k_over_h0 in np.logspace(-4.0, 6.0, 81):
                closure = curvature_constraint(
                    float(a),
                    float(k_over_h0),
                    params,
                    species,
                    xi_completion=xi_completion,
                )
                kinetic = np.asarray(closure["K"], dtype=float)
                mixing = np.asarray(closure["B"], dtype=float)
                mass = np.asarray(closure["M"], dtype=float)
                mixing_dot = np.asarray(closure["B_dot"], dtype=float)
                blocks = quadratic_blocks_with_time_derivative(
                    float(a),
                    float(k_over_h0),
                    params,
                    species,
                    xi_completion=xi_completion,
                )
                auxiliary = np.asarray(blocks["A"], dtype=float)
                kinetic_null = float(np.max(np.abs(kinetic[0, :])) / max(1.0, np.max(np.abs(kinetic))))
                antisymmetry = float(np.max(np.abs(mixing - mixing.T)) / max(1.0, np.max(np.abs(mixing))))
                mu_scale = max(1.0, abs(mass[0, 0]), abs(mixing_dot[0, 0]))
                mu_relative = abs(float(closure["mu_zeta"])) / mu_scale
                positive, negative, null, _ = normalized_inertia(kinetic[1:, 1:])
                maximum_kinetic_null = max(maximum_kinetic_null, kinetic_null)
                maximum_antisymmetry = max(maximum_antisymmetry, antisymmetry)
                minimum_q_tensor = min(minimum_q_tensor, float(blocks["rate"]["Q"]) + xi_completion)
                minimum_constraint_determinant = min(minimum_constraint_determinant, abs(float(np.linalg.det(auxiliary[:2, :2]))))
                minimum_mu_relative = min(minimum_mu_relative, mu_relative)
                rows.append(
                    {
                        "a": float(a),
                        "z": 1.0 / float(a) - 1.0,
                        "k_over_H0": float(k_over_h0),
                        "xi_completion": xi_completion,
                        "mu_zeta": float(closure["mu_zeta"]),
                        "mu_zeta_relative": mu_relative,
                        "constraint_determinant": float(np.linalg.det(auxiliary[:2, :2])),
                        "Q_tensor": float(blocks["rate"]["Q"]) + xi_completion,
                        "matter_kinetic_positive": positive,
                        "matter_kinetic_negative": negative,
                        "matter_kinetic_null": null,
                        "root_at_this_a": int(root is not None),
                        "root_k_over_H0": "" if root is None else root,
                    }
                )
        root_counts[xi_completion] = root_count

    output = ROOT / "generated" / "tables" / "rfg_xi_completion_audit.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if any(root_counts.values()):
        raise RuntimeError("Xi completion has a sampled mu_zeta root")
    if minimum_q_tensor <= 0.0 or minimum_constraint_determinant <= 0.0:
        raise RuntimeError("Xi completion has a tensor or auxiliary constraint singularity")
    if minimum_mu_relative <= 1.0e-8:
        raise RuntimeError("Xi completion approaches a sampled curvature-constraint singularity")
    if maximum_kinetic_null > 1.0e-10 or maximum_antisymmetry > 1.0e-8:
        raise RuntimeError("Xi completion breaks the exact DAE identities")

    print(f"Xi completion parameter                 = {XI_REFERENCE:.12g}")
    print(f"Xi completion cross-check               = {XI_CROSSCHECK:.12g}")
    print(f"sampled curvature-constraint roots      = {root_counts}")
    print(f"min Q_tensor                            = {minimum_q_tensor:.12e}")
    print(f"min |det A_(alpha,s)|                   = {minimum_constraint_determinant:.12e}")
    print(f"min normalized |mu_zeta|                = {minimum_mu_relative:.12e}")
    print(f"max transformed zeta kinetic residual   = {maximum_kinetic_null:.3e}")
    print(f"max reduced B antisymmetry residual      = {maximum_antisymmetry:.3e}")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
