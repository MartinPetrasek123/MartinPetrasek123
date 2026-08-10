#!/usr/bin/env python3
"""Action-level DAE closure audit for the RFG-R matter scalar sector.

The finite photon--baryon--CDM--massless-neutrino action has one exact null
kinetic direction.  This module identifies it before any hierarchy truncation:

    Delta_i = delta_i + 3 (1+w_i) zeta.

After eliminating the lapse, shift and fluid velocity auxiliaries, zeta has no
time derivative.  Its Euler equation is an algebraic curvature constraint.
All background-time derivatives used here are analytic action derivatives;
finite differences are used only by the separate validation script.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from extended_eft_mapping import extended_eft_w_coefficients
from multifluid_reduction import FluidSpecies, Planck2018Reference, quadratic_blocks
from rfg_regularized import RFGRegularizedParams


ROOT = Path(__file__).resolve().parents[1]


def density_curvature_transform(species: Iterable[FluidSpecies]) -> np.ndarray:
    """Map (zeta, Delta_i) to (zeta, delta_i) exactly."""
    sources = tuple(species)
    transform = np.eye(1 + len(sources), dtype=float)
    for index, source in enumerate(sources):
        transform[1 + index, 0] = -3.0 * (1.0 + source.w)
    return transform


def quadratic_blocks_with_time_derivative(
    a: float,
    k_over_h0: float,
    params: RFGRegularizedParams,
    species: Iterable[FluidSpecies],
    *,
    xi_completion: float = 0.0,
) -> dict[str, object]:
    """Return finite matrices and rates, optionally with Xi(R3+sigma^2).

    The Xi completion is a constant cosmological foliation coefficient.  Its
    operator vanishes exactly on FLRW, so it does not alter the reconstructed
    background.  The matrix increments below follow directly from the
    independently symbolically derived quadratic density.
    """
    sources = tuple(species)
    blocks = quadratic_blocks(a, k_over_h0, params, sources)
    rate = extended_eft_w_coefficients(a, params)
    h = float(rate["E"])
    k2 = k_over_h0 * k_over_h0
    gradient = k2 / (a * a)

    k0 = np.asarray(blocks["K0"], dtype=float)
    c = np.asarray(blocks["C"], dtype=float)
    auxiliary = np.asarray(blocks["A"], dtype=float)
    d = np.asarray(blocks["D"], dtype=float)
    m0 = np.asarray(blocks["M0"], dtype=float)
    k0_dot = np.zeros_like(k0)
    c_dot = np.zeros_like(c)
    auxiliary_dot = np.zeros_like(auxiliary)
    d_dot = np.zeros_like(d)
    m0_dot = np.zeros_like(m0)

    if not math.isfinite(xi_completion):
        raise ValueError("xi_completion must be finite")

    w0, w1 = float(rate["W0"]), float(rate["W1"])
    w4, w5, w6, w7 = (float(rate[key]) for key in ("W4", "W5", "W6", "W7"))
    w0_dot, w1_dot = float(rate["W0_dot"]), float(rate["W1_dot"])
    w4_dot, w5_dot, w6_dot, w7_dot = (float(rate[key]) for key in ("W4_dot", "W5_dot", "W6_dot", "W7_dot"))
    m5 = float(rate["m5_bar_hat"])
    m5_dot = float(rate["m5_bar_dot_hat"])

    k0_dot[0, 0] = -3.0 * a * a * (w5_dot + 2.0 * h * w5)
    c_dot[0, 0] = -3.0 * a * a * (w4_dot + 2.0 * h * w4)
    c_dot[1, 0] = -a * a * (w5_dot + 2.0 * h * w5)
    auxiliary_dot[0, 0] = 2.0 * w1_dot
    auxiliary_dot[0, 1] = auxiliary_dot[1, 0] = -a * a * (w4_dot + 2.0 * h * w4)
    auxiliary_dot[1, 1] = 2.0 * a**4 * (w7_dot + 4.0 * h * w7)
    d_dot[0, 0] = -w6_dot * k2
    d_dot[1, 0] = -2.0 * k2 * (m5_dot - 2.0 * h * m5) / (a * a)
    m0_dot[0, 0] = -2.0 * w0_dot * k2

    # For Xi(R3+sigma_ij sigma^ij), after the real-mode normalization used by
    # the finite action: Delta L/a^3 = Xi[2p alpha zeta + p zeta^2+s^2/3].
    # Xi is constant in the cosmological sector; the Weyl switch used for the
    # local GR domain is identically one on FLRW and is not differentiated here.
    xi_gradient = 2.0 * xi_completion * gradient
    auxiliary[1, 1] += 2.0 * xi_completion / 3.0
    d[0, 0] += xi_gradient
    m0[0, 0] += xi_gradient
    d_dot[0, 0] -= 2.0 * h * xi_gradient
    m0_dot[0, 0] -= 2.0 * h * xi_gradient

    for number, source in enumerate(sources):
        q_index, x_index = 1 + number, 2 + number
        rho = source.rho(a)
        rho_dot = -3.0 * h * (1.0 + source.w) * rho
        c_dot[x_index, 0] = -3.0 * rho_dot * (1.0 + source.w)
        c_dot[x_index, q_index] = -rho_dot
        auxiliary_dot[1, x_index] = auxiliary_dot[x_index, 1] = -rho_dot * (1.0 + source.w)
        auxiliary_dot[x_index, x_index] = -(1.0 + source.w) * gradient * (rho_dot - 2.0 * h * rho)
        d_dot[0, q_index] = -rho_dot
        m0_dot[q_index, q_index] = -rho_dot * source.w / (1.0 + source.w)

    return {
        **blocks,
        "K0_dot": k0_dot,
        "C_dot": c_dot,
        "A_dot": auxiliary_dot,
        "D_dot": d_dot,
        "M0_dot": m0_dot,
        "rate": rate,
        "xi_completion": xi_completion,
    }


def reduced_density_blocks(
    a: float,
    k_over_h0: float,
    params: RFGRegularizedParams,
    species: Iterable[FluidSpecies],
    *,
    xi_completion: float = 0.0,
) -> dict[str, np.ndarray | float]:
    """Eliminate auxiliaries and transform to (zeta, Delta_i)."""
    sources = tuple(species)
    blocks = quadratic_blocks_with_time_derivative(
        a,
        k_over_h0,
        params,
        sources,
        xi_completion=xi_completion,
    )
    k0, c, auxiliary, d, m0 = (np.asarray(blocks[key], dtype=float) for key in ("K0", "C", "A", "D", "M0"))
    k0_dot, c_dot, auxiliary_dot, d_dot, m0_dot = (
        np.asarray(blocks[key], dtype=float) for key in ("K0_dot", "C_dot", "A_dot", "D_dot", "M0_dot")
    )
    inverse = np.linalg.inv(auxiliary)
    kinetic = k0 - c.T @ inverse @ c
    mixing = -c.T @ inverse @ d
    mass = m0 - d.T @ inverse @ d
    kinetic_dot = k0_dot - c_dot.T @ inverse @ c - c.T @ inverse @ c_dot + c.T @ inverse @ auxiliary_dot @ inverse @ c
    mixing_dot = -c_dot.T @ inverse @ d + c.T @ inverse @ auxiliary_dot @ inverse @ d - c.T @ inverse @ d_dot
    mass_dot = m0_dot - d_dot.T @ inverse @ d - d.T @ inverse @ d_dot + d.T @ inverse @ auxiliary_dot @ inverse @ d
    transform = density_curvature_transform(sources)
    return {
        "K": transform.T @ kinetic @ transform,
        "B": transform.T @ mixing @ transform,
        "M": transform.T @ mass @ transform,
        "K_dot": transform.T @ kinetic_dot @ transform,
        "B_dot": transform.T @ mixing_dot @ transform,
        "M_dot": transform.T @ mass_dot @ transform,
        "transform": transform,
        "hubble": float(blocks["rate"]["E"]),
    }


def curvature_constraint(
    a: float,
    k_over_h0: float,
    params: RFGRegularizedParams,
    species: Iterable[FluidSpecies],
    *,
    xi_completion: float = 0.0,
) -> dict[str, np.ndarray | float]:
    """Return the zeta Euler constraint after integrating its dot(zeta) term.

    For L=dot(y)^T K dot(y)/2+dot(y)^T B y+y^T M y/2 and
    y=(zeta,Delta_i), K_{0i}=0.  The zeta equation is

      mu_zeta zeta + c_i Delta_i + u_i dot(Delta_i)=0,

    where the final velocity coefficient is retained even though it vanishes
    identically for this action (B is symmetric).
    """
    reduced = reduced_density_blocks(
        a,
        k_over_h0,
        params,
        species,
        xi_completion=xi_completion,
    )
    mixing = np.asarray(reduced["B"], dtype=float)
    mass = np.asarray(reduced["M"], dtype=float)
    mixing_dot = np.asarray(reduced["B_dot"], dtype=float)
    return {
        **reduced,
        "mu_zeta": float(mass[0, 0] - mixing_dot[0, 0]),
        "density_coefficient": mass[0, 1:] - mixing_dot[0, 1:],
        "velocity_coefficient": mixing[1:, 0] - mixing[0, 1:],
    }


def _relative_scale(values: np.ndarray) -> float:
    return float(np.max(np.abs(values))) if values.size else 0.0


def find_mu_root(
    a: float,
    params: RFGRegularizedParams,
    species: Iterable[FluidSpecies],
    *,
    k_min: float = 1.0e-4,
    k_max: float = 1.0e6,
    samples: int = 401,
    xi_completion: float = 0.0,
) -> float | None:
    """Locate a sign-changing analytic curvature-constraint coefficient."""
    grid = np.logspace(math.log10(k_min), math.log10(k_max), samples)
    left = float(grid[0])
    left_value = float(curvature_constraint(a, left, params, species, xi_completion=xi_completion)["mu_zeta"])
    for current in grid[1:]:
        right = float(current)
        right_value = float(curvature_constraint(a, right, params, species, xi_completion=xi_completion)["mu_zeta"])
        if left_value == 0.0:
            return left
        if left_value * right_value < 0.0:
            for _ in range(96):
                middle = math.sqrt(left * right)
                middle_value = float(curvature_constraint(a, middle, params, species, xi_completion=xi_completion)["mu_zeta"])
                if left_value * middle_value <= 0.0:
                    right, right_value = middle, middle_value
                else:
                    left, left_value = middle, middle_value
            return math.sqrt(left * right)
        left, left_value = right, right_value
    return None


def main() -> None:
    reference = Planck2018Reference()
    params = reference.rfg_background_parameters()
    sources = reference.species()
    rows: list[dict[str, float | str]] = []
    maximum_kinetic_null = 0.0
    maximum_antisymmetry = 0.0
    for a in np.logspace(-7.0, 0.0, 49):
        reduced = curvature_constraint(float(a), 1.0, params, sources)
        kinetic = np.asarray(reduced["K"], dtype=float)
        mixing = np.asarray(reduced["B"], dtype=float)
        kinetic_null = _relative_scale(kinetic[0, :]) / max(1.0, _relative_scale(kinetic))
        antisymmetry = _relative_scale(mixing - mixing.T) / max(1.0, _relative_scale(mixing))
        maximum_kinetic_null = max(maximum_kinetic_null, kinetic_null)
        maximum_antisymmetry = max(maximum_antisymmetry, antisymmetry)
        root = find_mu_root(float(a), params, sources)
        rows.append(
            {
                "a": float(a),
                "z": 1.0 / float(a) - 1.0,
                "kinetic_zeta_relative_residual": kinetic_null,
                "mixing_antisymmetry_relative_residual": antisymmetry,
                "mu_zeta_at_k_over_H0_1": float(reduced["mu_zeta"]),
                "mu_zeta_root_k_over_H0": "" if root is None else root,
            }
        )

    output = ROOT / "generated" / "tables" / "multifluid_dae_closure.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    roots = [row["mu_zeta_root_k_over_H0"] for row in rows if row["mu_zeta_root_k_over_H0"] != ""]
    print(f"max transformed zeta kinetic residual = {maximum_kinetic_null:.3e}")
    print(f"max reduced B antisymmetry residual    = {maximum_antisymmetry:.3e}")
    print(f"curvature-constraint roots on grid     = {len(roots)} / {len(rows)}")
    if roots:
        print(f"latest root k/H0                       = {float(roots[-1]):.12g}")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
