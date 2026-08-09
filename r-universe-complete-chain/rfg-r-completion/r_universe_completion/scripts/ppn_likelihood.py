#!/usr/bin/env python3
"""Local GR-matching and the Gaussian Cassini gamma likelihood."""

from __future__ import annotations

import math


# Dimensionless matching thresholds W = sqrt(C_abcd C^abcd)/(H0/c)^2.
W_COSMOLOGY = 1.0e8
W_LOCAL_GR = 1.0e9


def bump(x: float) -> float:
    return 0.0 if x <= 0.0 else math.exp(-1.0 / x)


def cosmology_weight(weyl_ratio: float) -> float:
    """C-infinity interpolation: one in the cosmological EFT, zero in local GR."""
    if weyl_ratio <= W_COSMOLOGY:
        return 1.0
    if weyl_ratio >= W_LOCAL_GR:
        return 0.0
    t = (weyl_ratio - W_COSMOLOGY) / (W_LOCAL_GR - W_COSMOLOGY)
    left = bump(1.0 - t)
    right = bump(t)
    return left / (left + right)


def solar_weyl_ratio(radius_m: float, h0_km_s_mpc: float = 67.4) -> float:
    """W for a Schwarzschild Sun, using sqrt(C^2)=sqrt(48) GM/(c^2 r^3)."""
    gm_over_c2_m = 1476.6250385
    megaparsec_m = 3.085677581491367e22
    c_m_s = 299792458.0
    h0_s = h0_km_s_mpc * 1000.0 / megaparsec_m
    curvature = math.sqrt(48.0) * gm_over_c2_m / radius_m**3
    h0_curvature = (h0_s / c_m_s) ** 2
    return curvature / h0_curvature


def ppn_parameters(weyl_ratio: float) -> dict[str, float]:
    """The matched local domain has the Einstein-Hilbert action exactly."""
    if cosmology_weight(weyl_ratio) != 0.0:
        raise ValueError("PPN prediction is defined only inside the local GR-matched domain")
    return {"gamma": 1.0, "beta": 1.0, "alpha1": 0.0, "alpha2": 0.0}


def cassini_minus2loglike(gamma: float) -> float:
    """Bertotti, Iess and Tortora (2003): gamma-1=(2.1 +- 2.3)e-5."""
    return ((gamma - 1.0 - 2.1e-5) / 2.3e-5) ** 2


def main() -> None:
    au_m = 149597870700.0
    ratio = solar_weyl_ratio(au_m)
    params = ppn_parameters(ratio)
    print(f"Solar Weyl ratio at 1 AU = {ratio:.6e}")
    print(f"Cosmological action weight = {cosmology_weight(ratio):.1f}")
    print("PPN: gamma={gamma:.1f}, beta={beta:.1f}, alpha1={alpha1:.1f}, alpha2={alpha2:.1f}".format(**params))
    print(f"Cassini -2 ln L = {cassini_minus2loglike(params['gamma']):.10f}")


if __name__ == "__main__":
    main()

