#!/usr/bin/env python3
"""Independent checks for the exact RFG-R multi-species reduction."""

from __future__ import annotations

import csv
from fractions import Fraction
import math
from pathlib import Path

import numpy as np

from multifluid_reduction import (
    FluidSpecies,
    Planck2018Reference,
    ScalarMetricSource,
    baryon_cdm_continuity_rhs,
    massless_hierarchy_rhs,
    massless_anisotropic_stress,
    normalized_inertia,
    photon_collision,
    quadratic_blocks,
    reduce_auxiliaries,
    spatial_traceless_shear_rhs_from_coefficients,
)


ROOT = Path(__file__).resolve().parents[1]


def _fraction_zero(rows: int, columns: int) -> list[list[Fraction]]:
    return [[Fraction(0) for _ in range(columns)] for _ in range(rows)]


def _fraction_inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    augmented = [row[:] + [Fraction(int(i == j)) for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row][column] != 0), None)
        if pivot is None:
            raise ZeroDivisionError("singular exact matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [entry / scale for entry in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [entry - factor * reference for entry, reference in zip(augmented[row], augmented[column])]
    return [row[size:] for row in augmented]


def _fraction_multiply(left: list[list[Fraction]], right: list[list[Fraction]]) -> list[list[Fraction]]:
    return [[sum(left[row][index] * right[index][column] for index in range(len(right))) for column in range(len(right[0]))] for row in range(len(left))]


def _fraction_transpose(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    return [list(column) for column in zip(*matrix)]


def _fraction_rank(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    rank = 0
    columns = len(work[0])
    for column in range(columns):
        pivot = next((row for row in range(rank, len(work)) if work[row][column] != 0), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][column]
        work[rank] = [entry / scale for entry in work[rank]]
        for row in range(len(work)):
            if row != rank and work[row][column] != 0:
                factor = work[row][column]
                work[row] = [entry - factor * reference for entry, reference in zip(work[row], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def _exact_gr_kinetic_rank() -> int:
    """Do the GR check with rational arithmetic: no tolerance enters the rank."""
    # a=H=k=1 and rho_b,rho_c,rho_gamma,rho_nu=(3,6,9,12)/10.  The sum is
    # 3H^2 exactly.  The W_i are the exact GR values in the conventions used
    # by the RFG-R extended action.
    densities = (Fraction(3, 10), Fraction(6, 10), Fraction(9, 10), Fraction(12, 10))
    equations_of_state = (Fraction(0), Fraction(0), Fraction(1, 3), Fraction(1, 3))
    sources = len(densities)
    n, m = 1 + sources, 2 + sources
    k0, c, a_matrix, _, _ = (_fraction_zero(n, n), _fraction_zero(m, n), _fraction_zero(m, m), _fraction_zero(m, n), _fraction_zero(n, n))
    # W1=-3, W4=-2, W5=2, W7=0.  The remaining entries do not enter K.
    k0[0][0] = Fraction(-6)
    c[0][0], c[1][0] = Fraction(6), Fraction(-2)
    a_matrix[0][0], a_matrix[0][1], a_matrix[1][0], a_matrix[1][1] = Fraction(-6), Fraction(2), Fraction(2), Fraction(0)
    for index, (rho, w) in enumerate(zip(densities, equations_of_state)):
        q_index, x_index = 1 + index, 2 + index
        c[x_index][0] = -3 * rho * (1 + w)
        c[x_index][q_index] = -rho
        a_matrix[1][x_index] = a_matrix[x_index][1] = -rho * (1 + w)
        a_matrix[x_index][x_index] = -rho * (1 + w)

    inverse = _fraction_inverse(a_matrix)
    correction = _fraction_multiply(_fraction_multiply(_fraction_transpose(c), inverse), c)
    kinetic = [[k0[row][column] - correction[row][column] for column in range(n)] for row in range(n)]
    return _fraction_rank(kinetic)


def _check_hierarchy_action_match() -> None:
    """The photon ell=0 equation must reproduce the action continuity equation."""
    metric = ScalarMetricSource(alpha=0.17, zeta_dot=-0.23, s=0.41)
    k_over_a = 1.9
    theta_1 = -0.37
    theta_0_dot = massless_hierarchy_rhs(0, 0.0, 0.0, theta_1, metric, k_over_a)
    delta_dot_from_hierarchy = 4.0 * theta_0_dot
    delta_dot_from_conservation = -4.0 * k_over_a * theta_1 / 3.0 - 4.0 * metric.zeta_dot - 4.0 * metric.s / 3.0
    assert abs(delta_dot_from_hierarchy - delta_dot_from_conservation) < 1.0e-14
    assert photon_collision(0, 0.2, 0.3, -0.1, 0.4, 7.0) == 0.0
    assert abs(baryon_cdm_continuity_rhs(metric.zeta_dot, metric.s, theta_1, k_over_a) - (-k_over_a * theta_1 - 3.0 * metric.zeta_dot - metric.s)) < 1.0e-14


def _check_spatial_traceless_equation() -> None:
    """Verify the new shear equation exactly reduces to its GR identity."""
    a, k, h = 0.73, 1.9, 0.81
    alpha, zeta, zeta_dot, shear, pi = 0.11, -0.23, 0.37, -0.19, 0.07
    actual = spatial_traceless_shear_rhs_from_coefficients(
        a=a,
        k_over_h0=k,
        hubble=h,
        q=1.0,
        q_x=0.0,
        q_dot=0.0,
        alpha=alpha,
        zeta=zeta,
        zeta_dot=zeta_dot,
        shear=shear,
        anisotropic_stress=pi,
    )
    expected = -3.0 * h * shear - (k * k / (a * a)) * (alpha + zeta) - pi
    assert abs(actual - expected) < 1.0e-14
    assert massless_anisotropic_stress(3.0, 0.25) == 0.6


def main() -> None:
    assert _exact_gr_kinetic_rank() == 4
    _check_hierarchy_action_match()
    _check_spatial_traceless_equation()

    reference = Planck2018Reference()
    params = reference.rfg_background_parameters()
    params.validate()
    species = reference.species()
    assert {source.name for source in species} == {
        "baryon",
        "cdm",
        "photon_monopole_dipole",
        "massless_neutrino_monopole_dipole",
    }
    assert abs(params.omega_m0 - 0.3137721026737606) < 1.0e-15
    assert abs(params.omega_r0 - 9.219892755013364e-5) < 1.0e-19

    rows: list[dict[str, float | int]] = []
    maximum_constraint_error = 0.0
    maximum_null_residual = 0.0
    for scale_factor in np.logspace(-7.0, 0.0, 25):
        for k_over_h0 in np.logspace(-2.0, 2.0, 17):
            blocks = quadratic_blocks(float(scale_factor), float(k_over_h0), params, species)
            reduced = reduce_auxiliaries(blocks)
            kinetic = reduced["K"]
            assert np.allclose(kinetic, kinetic.T, rtol=0.0, atol=2.0e-12 * max(1.0, float(np.max(np.abs(kinetic)))))
            assert blocks["constraint_discriminant"] > 0.0
            exact_constraint_determinant = -(scale_factor**4) * blocks["constraint_discriminant"]
            constraint_error = abs(blocks["constraint_block_determinant"] - exact_constraint_determinant) / max(1.0, abs(exact_constraint_determinant))
            maximum_constraint_error = max(maximum_constraint_error, constraint_error)
            positive, negative, null, null_residual = normalized_inertia(kinetic)
            maximum_null_residual = max(maximum_null_residual, null_residual)
            assert (positive, negative, null) == (4, 0, 1)
            rows.append(
                {
                    "a": float(scale_factor),
                    "z": 1.0 / float(scale_factor) - 1.0,
                    "k_over_H0": float(k_over_h0),
                    "constraint_discriminant": float(blocks["constraint_discriminant"]),
                    "constraint_block_determinant": float(blocks["constraint_block_determinant"]),
                    "constraint_relative_error": constraint_error,
                    "kinetic_positive": positive,
                    "kinetic_negative": negative,
                    "kinetic_null": null,
                    "normalized_null_residual": null_residual,
                    "m5_bar_hat": float(blocks["m5_bar_hat"]),
                }
            )

    assert maximum_constraint_error < 3.0e-14
    assert maximum_null_residual < 3.0e-10
    assert max(abs(float(row["m5_bar_hat"])) for row in rows) > 1.0e-6
    output = ROOT / "generated" / "tables" / "multifluid_core_audit.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print("Exact GR rational multi-fluid rank     = 4")
    print("GR physical content                    = four matter scalar modes, no gravity scalar")
    print(f"Planck reference grid points           = {len(rows)}")
    print(f"max lapse-shift determinant residual   = {maximum_constraint_error:.3e}")
    print(f"max normalized kinetic null residual   = {maximum_null_residual:.3e}")
    print("RFG-R core inertia                     = (positive, negative, null) = (4, 0, 1)")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
