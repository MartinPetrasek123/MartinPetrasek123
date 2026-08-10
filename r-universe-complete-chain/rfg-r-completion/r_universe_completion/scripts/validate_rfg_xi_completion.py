#!/usr/bin/env python3
"""Regression tests for the Xi(R3+sigma^2) action completion."""

from __future__ import annotations

import sympy as sp

from derive_spatial_traceless_equation import quadratic_mode_lagrangian
from multifluid_reduction import Planck2018Reference
from rfg_dae_closure import curvature_constraint, find_mu_root
from rfg_xi_completion import XI_AUDIT_VALUES, XI_REFERENCE


def main() -> None:
    value, lagrangian = quadratic_mode_lagrangian()
    fixed = {
        value.shear: 0,
        value.shear_dot: 0,
        value.eta: 0,
        value.sigma_completion: 0,
        value.geometry_completion: 0,
    }
    base = sp.factor(lagrangian.subs(fixed))
    completed = sp.factor(
        lagrangian.subs(
            {
                value.shear: 0,
                value.shear_dot: 0,
                value.eta: 0,
                value.sigma_completion: 0,
            }
        )
    )
    expected = value.geometry_completion * value.k**2 * (
        6 * value.a**2 * value.alpha * value.zeta
        + 3 * value.a**2 * value.zeta**2
        + value.beta**2 * value.k**2
    ) / (3 * value.a)
    assert sp.simplify(completed - base - expected) == 0

    reference = Planck2018Reference()
    params = reference.rfg_background_parameters()
    species = reference.species()
    assert find_mu_root(1.0, params, species) is not None
    for xi_completion in XI_AUDIT_VALUES:
        assert find_mu_root(1.0, params, species, xi_completion=xi_completion) is None
        closure = curvature_constraint(1.0, 1.0, params, species, xi_completion=xi_completion)
        assert float(closure["mu_zeta"]) > 0.0

    print("Xi completion symbolic and DAE regression checks passed")
    print("exact Xi quadratic density recovered")
    print(f"a=1 Xi={XI_REFERENCE:g} and Xi={XI_AUDIT_VALUES[-1]:g} have no sampled mu_zeta root")


if __name__ == "__main__":
    main()
