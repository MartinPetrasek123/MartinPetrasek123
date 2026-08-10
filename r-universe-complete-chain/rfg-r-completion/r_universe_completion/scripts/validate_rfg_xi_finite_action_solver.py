#!/usr/bin/env python3
"""Regression checks for the implicit finite-action RFG-RXi DAE solver."""

import numpy as np

from extended_eft_mapping import extended_eft_coefficients
from multifluid_reduction import FluidSpecies
from rfg_xi_finite_action_solver import GRLimitParams, integrate_finite_action


def main() -> None:
    solutions = {
        xi_completion: integrate_finite_action(k_over_h0=1.0, xi_completion=xi_completion, samples=41)
        for xi_completion in (1.0, 2.0)
    }
    for solution in solutions.values():
        assert np.all(np.isfinite(solution.state))
        assert np.all(np.isfinite(solution.zeta))
        assert np.max(solution.constraint_relative_residual) < 1.0e-16
        assert solution.minimum_mu_relative > 1.0e-8

    unit = integrate_finite_action(k_over_h0=1.0, xi_completion=1.0, samples=41, seed_amplitude=1.0e-8)
    doubled = integrate_finite_action(k_over_h0=1.0, xi_completion=1.0, samples=41, seed_amplitude=2.0e-8)
    state_scale_error = np.max(np.abs(doubled.state - 2.0 * unit.state)) / np.max(np.abs(doubled.state))
    zeta_scale_error = np.max(np.abs(doubled.zeta - 2.0 * unit.zeta)) / np.max(np.abs(doubled.zeta))
    assert state_scale_error < 1.0e-10
    assert zeta_scale_error < 1.0e-10

    gr_params = GRLimitParams()
    for a in (1.0e-8, 1.0e-4, 1.0):
        coefficients = extended_eft_coefficients(a, gr_params)
        assert coefficients["Q"] == 1.0
        assert coefficients["Q_X"] == 0.0
        assert coefficients["m5_bar_hat"] == 0.0
    gr_species = (
        FluidSpecies("baryon", 0.15, 0.0),
        FluidSpecies("cdm", 0.8499, 0.0),
        FluidSpecies("photon", 0.00005, 1.0 / 3.0),
        FluidSpecies("massless_neutrino", 0.00005, 1.0 / 3.0),
    )
    gr_solution = integrate_finite_action(
        k_over_h0=1.0,
        xi_completion=0.0,
        samples=41,
        params=gr_params,
        species=gr_species,
    )
    assert np.all(np.isfinite(gr_solution.state))
    assert np.max(gr_solution.constraint_relative_residual) < 1.0e-16
    print("RFG-RXi finite-action implicit DAE solver checks passed")
    print(f"max constraint residual = {max(np.max(solution.constraint_relative_residual) for solution in solutions.values()):.3e}")
    print(f"linearity scale errors  = state {state_scale_error:.3e}, zeta {zeta_scale_error:.3e}")
    print("exact GR action limit   = Q=1, Q_X=m5_bar=0; finite DAE trajectory passed")


if __name__ == "__main__":
    main()
