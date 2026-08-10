#!/usr/bin/env python3
"""Implicitly integrate the finite action-level RFG-RXi scalar DAE reduction.

This is a solver kernel, not a CMB code: photons and neutrinos are retained
only through their monopole/dipole variables already present in the exact
finite action.  The untruncated kinetic hierarchies and visibility/line-of-
sight calculation are a subsequent extension.  No Einstein equation or
quasi-static effective coupling is substituted for the RFG-RXi constraints.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import partial
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.integrate import solve_ivp

from extended_eft_mapping import background_derivatives
from multifluid_reduction import FluidSpecies, Planck2018Reference
from rfg_dae_closure import curvature_constraint
from rfg_regularized import RFGRegularizedParams


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "generated" / "tables"


@dataclass(frozen=True)
class FiniteActionSolution:
    """A converged finite-action DAE trajectory in log-scale-factor time."""

    xi_completion: float
    k_over_h0: float
    a: np.ndarray
    state: np.ndarray
    zeta: np.ndarray
    constraint_relative_residual: np.ndarray
    minimum_mu_relative: float


@dataclass(frozen=True)
class GRLimitParams:
    """Duck-typed exact Einstein-Hilbert limit of the regular action."""

    omega_m0: float = 0.9999
    omega_r0: float = 0.0001
    theta: float = 1.6
    epsilon: float = 1.0e-8
    p: int = 4

    @property
    def omega_R0(self) -> float:
        return 0.0

    @property
    def A(self) -> float:
        return 0.0

    def validate(self) -> None:
        if not (0.0 < self.theta < 2.0) or self.p <= self.theta or self.p % 2 != 0:
            raise ValueError("invalid GR-limit regularization parameters")


def _constraint_zeta(closure: dict[str, object], delta: np.ndarray) -> tuple[float, float]:
    """Solve the exact algebraic curvature constraint for zeta."""
    mu = float(closure["mu_zeta"])
    density = np.asarray(closure["density_coefficient"], dtype=float)
    scale = max(1.0, abs(mu), float(np.max(np.abs(density))))
    if abs(mu) / scale <= 1.0e-10:
        raise RuntimeError("RFG-RXi curvature constraint is singular on the integration path")
    zeta = -float(density @ delta) / mu
    residual = abs(mu * zeta + float(density @ delta)) / scale
    return zeta, residual


def finite_action_rhs(
    log_a: float,
    state: np.ndarray,
    *,
    k_over_h0: float,
    params: RFGRegularizedParams,
    species: Iterable[FluidSpecies],
    xi_completion: float,
) -> np.ndarray:
    """Return d(Delta_i, dDelta_i/dln a)/dln a from the reduced action.

    For y=(zeta,Delta_i), the exact reduced density has K_0i=0 and symmetric
    B.  The zeta row is therefore algebraic.  The material rows obey

      K_ij Delta_j,tt + Kdot_ij Delta_j,t + (Bdot-M)_iA y_A = 0,

    with all coefficients computed as analytic derivatives of the action.
    """
    a = math.exp(log_a)
    sources = tuple(species)
    count = len(sources)
    delta = np.asarray(state[:count], dtype=float)
    delta_log_derivative = np.asarray(state[count:], dtype=float)
    closure = curvature_constraint(a, k_over_h0, params, sources, xi_completion=xi_completion)
    zeta, _ = _constraint_zeta(closure, delta)
    kinetic = np.asarray(closure["K"], dtype=float)[1:, 1:]
    kinetic_dot = np.asarray(closure["K_dot"], dtype=float)[1:, 1:]
    force = np.asarray(closure["B_dot"], dtype=float) - np.asarray(closure["M"], dtype=float)
    hubble = float(closure["hubble"])
    h_dot = float(background_derivatives(a, params)["Hdot_over_H0_sq"])
    material_force = force[1:, 0] * zeta + force[1:, 1:] @ delta
    acceleration = np.linalg.solve(
        kinetic,
        -((h_dot / (hubble * hubble)) * kinetic + kinetic_dot / hubble) @ delta_log_derivative
        - material_force / (hubble * hubble),
    )
    return np.concatenate((delta_log_derivative, acceleration))


def integrate_finite_action(
    *,
    k_over_h0: float,
    xi_completion: float,
    a_start: float = 1.0e-8,
    a_end: float = 1.0,
    samples: int = 81,
    seed_amplitude: float = 1.0e-8,
    params: RFGRegularizedParams | None = None,
    species: Iterable[FluidSpecies] | None = None,
) -> FiniteActionSolution:
    """Integrate one finite-action RFG-RXi DAE trajectory with Radau."""
    if not (0.0 < a_start < a_end <= 1.0) or samples < 2:
        raise ValueError("require 0 < a_start < a_end <= 1 and at least two samples")
    reference = Planck2018Reference()
    selected_params = params if params is not None else reference.rfg_background_parameters()
    sources = tuple(species) if species is not None else reference.species()
    count = len(sources)
    # A nonzero linear seed exercises the DAE; it is deliberately not yet
    # labelled an adiabatic primordial mode.
    seed = seed_amplitude * np.array((3.0, 3.0, 4.0, 4.0), dtype=float)
    if count != len(seed):
        raise ValueError("the finite action solver currently expects four material species")
    initial = np.concatenate((seed, np.zeros(count, dtype=float)))
    log_grid = np.linspace(math.log(a_start), math.log(a_end), samples)
    rhs = partial(
        finite_action_rhs,
        k_over_h0=k_over_h0,
        params=selected_params,
        species=sources,
        xi_completion=xi_completion,
    )
    result = solve_ivp(
        rhs,
        (float(log_grid[0]), float(log_grid[-1])),
        initial,
        method="Radau",
        t_eval=log_grid,
        rtol=2.0e-8,
        atol=seed_amplitude * 1.0e-5,
        first_step=1.0e-4,
        max_step=5.0e-2,
    )
    if not result.success:
        raise RuntimeError(f"RFG-RXi finite-action integration failed: {result.message}")
    zeta_values: list[float] = []
    residual_values: list[float] = []
    mu_relative_values: list[float] = []
    for log_a, column in zip(result.t, result.y.T):
        closure = curvature_constraint(
            math.exp(float(log_a)), k_over_h0, selected_params, sources, xi_completion=xi_completion
        )
        zeta, residual = _constraint_zeta(closure, column[:count])
        density = np.asarray(closure["density_coefficient"], dtype=float)
        mu_relative_values.append(abs(float(closure["mu_zeta"])) / max(1.0, abs(float(closure["mu_zeta"])), float(np.max(np.abs(density)))))
        zeta_values.append(zeta)
        residual_values.append(residual)
    return FiniteActionSolution(
        xi_completion=xi_completion,
        k_over_h0=k_over_h0,
        a=np.exp(result.t),
        state=result.y.T,
        zeta=np.asarray(zeta_values),
        constraint_relative_residual=np.asarray(residual_values),
        minimum_mu_relative=float(min(mu_relative_values)),
    )


def write_solutions(solutions: dict[float, FiniteActionSolution]) -> Path:
    """Write finite-action DAE trajectories without labelling them spectra."""
    rows: list[dict[str, float]] = []
    for xi_completion, solution in solutions.items():
        for index, a in enumerate(solution.a):
            row: dict[str, float] = {
                "xi_completion": xi_completion,
                "k_over_H0": solution.k_over_h0,
                "a": float(a),
                "z": 1.0 / float(a) - 1.0,
                "zeta": float(solution.zeta[index]),
                "constraint_relative_residual": float(solution.constraint_relative_residual[index]),
            }
            for species_index in range(solution.state.shape[1] // 2):
                row[f"Delta_{species_index}"] = float(solution.state[index, species_index])
                row[f"dDelta_dln_a_{species_index}"] = float(solution.state[index, species_index + solution.state.shape[1] // 2])
            rows.append(row)
    output = TABLES / "rfg_xi_finite_action_solver.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return output


def main() -> None:
    solutions: dict[float, FiniteActionSolution] = {}
    for xi_completion in (1.0, 2.0):
        solution = integrate_finite_action(k_over_h0=1.0, xi_completion=xi_completion)
        solutions[xi_completion] = solution
        print(
            f"Xi={xi_completion:g}: points={len(solution.a)}, "
            f"max constraint residual={np.max(solution.constraint_relative_residual):.3e}, "
            f"min normalized |mu_zeta|={solution.minimum_mu_relative:.3e}"
        )
    print(f"wrote {write_solutions(solutions)}")


if __name__ == "__main__":
    main()
