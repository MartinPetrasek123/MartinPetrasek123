#!/usr/bin/env python3
"""Independent numerical checks of the analytic RFG-R DAE reduction."""

from __future__ import annotations

import math

import numpy as np

from multifluid_reduction import Planck2018Reference
from rfg_dae_closure import curvature_constraint, find_mu_root, reduced_density_blocks


def main() -> None:
    reference = Planck2018Reference()
    params = reference.rfg_background_parameters()
    species = reference.species()
    maximum_kinetic_null = 0.0
    maximum_antisymmetry = 0.0
    maximum_rate_error = 0.0
    step = 1.0e-5

    for a, k in ((1.0e-4, 1.0e-2), (1.0e-2, 1.0), (0.1, 10.0), (0.5, 1.0), (1.0, 1.0)):
        reduced = curvature_constraint(a, k, params, species)
        kinetic = np.asarray(reduced["K"], dtype=float)
        mixing = np.asarray(reduced["B"], dtype=float)
        kinetic_null = np.max(np.abs(kinetic[0, :])) / max(1.0, np.max(np.abs(kinetic)))
        antisymmetry = np.max(np.abs(mixing - mixing.T)) / max(1.0, np.max(np.abs(mixing)))
        maximum_kinetic_null = max(maximum_kinetic_null, float(kinetic_null))
        maximum_antisymmetry = max(maximum_antisymmetry, float(antisymmetry))

        plus = reduced_density_blocks(a * math.exp(step), k, params, species)
        minus = reduced_density_blocks(a * math.exp(-step), k, params, species)
        h = float(reduced["hubble"])
        numerical_rate = h * (np.asarray(plus["B"]) - np.asarray(minus["B"])) / (2.0 * step)
        analytic_rate = np.asarray(reduced["B_dot"], dtype=float)
        rate_error = np.max(np.abs(numerical_rate - analytic_rate)) / max(1.0, np.max(np.abs(numerical_rate)), np.max(np.abs(analytic_rate)))
        maximum_rate_error = max(maximum_rate_error, float(rate_error))

    root = find_mu_root(1.0, params, species)
    assert root is not None
    at_root = curvature_constraint(1.0, root, params, species)
    left = float(curvature_constraint(1.0, root * (1.0 - 1.0e-5), params, species)["mu_zeta"])
    right = float(curvature_constraint(1.0, root * (1.0 + 1.0e-5), params, species)["mu_zeta"])
    root_residual = abs(float(at_root["mu_zeta"])) / max(1.0, abs(left), abs(right))

    assert maximum_kinetic_null < 1.0e-11
    assert maximum_antisymmetry < 1.0e-11
    assert maximum_rate_error < 2.0e-7
    assert left * right < 0.0
    assert root_residual < 1.0e-8

    print("RFG-R DAE closure derivation validation OK")
    print(f"max transformed zeta kinetic residual = {maximum_kinetic_null:.3e}")
    print(f"max B antisymmetry residual            = {maximum_antisymmetry:.3e}")
    print(f"max analytic B-dot check residual      = {maximum_rate_error:.3e}")
    print(f"a=1 curvature-constraint root k/H0     = {root:.12g}")
    print(f"normalized root residual                = {root_residual:.3e}")


if __name__ == "__main__":
    main()
