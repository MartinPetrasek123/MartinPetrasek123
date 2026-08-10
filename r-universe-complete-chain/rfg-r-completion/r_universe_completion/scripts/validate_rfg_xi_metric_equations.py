#!/usr/bin/env python3
"""Regress the full scalar RFG-RXi metric variation before any solver use."""

import math

import sympy as sp

from derive_spatial_traceless_equation import spatial_gauge_metric_residuals
from multifluid_reduction import spatial_traceless_shear_rhs_from_coefficients


def main() -> None:
    value, residuals = spatial_gauge_metric_residuals()
    base = {value.geometry_completion: 0, value.sigma_completion: 0}
    completed = {value.sigma_completion: 0}
    expected = {
        "lapse": 2 * value.geometry_completion * value.a * value.k**2 * value.zeta,
        "shift": 2 * value.geometry_completion * value.beta * value.k**4 / (3 * value.a),
        "trace": -2 * value.geometry_completion * value.a * value.k**2 * (value.alpha + value.zeta),
        "traceless": -2
        * value.geometry_completion
        * value.a
        * value.k**4
        * (value.h * value.beta + value.alpha + value.beta_dot + value.zeta)
        / 3,
    }
    for name, identity in expected.items():
        delta = sp.factor(residuals[name].subs(completed) - residuals[name].subs(base))
        assert sp.simplify(delta - identity) == 0, name

    common = {
        "a": 0.37,
        "k_over_h0": 0.81,
        "hubble": 2.4,
        "q_x": -0.13,
        "q_dot": 0.27,
        "alpha": 0.02,
        "zeta": -0.03,
        "zeta_dot": 0.04,
        "shear": -0.05,
        "anisotropic_stress": 0.06,
    }
    completed_value = spatial_traceless_shear_rhs_from_coefficients(
        q=0.73, xi_completion=1.0, **common
    )
    direct_shift_value = spatial_traceless_shear_rhs_from_coefficients(
        q=1.73, xi_completion=0.0, **common
    )
    assert math.isclose(completed_value, direct_shift_value, rel_tol=0.0, abs_tol=1.0e-14)
    print("RFG-RXi lapse, shift, trace, and traceless variations reproduced")
    print("RFG-RXi shear equation uses the exact Q -> Q+Xi completion")


if __name__ == "__main__":
    main()
