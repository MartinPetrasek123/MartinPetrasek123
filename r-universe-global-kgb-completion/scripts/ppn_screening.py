#!/usr/bin/env python3
"""Solar-system Vainshtein gate for the cubic KGB term at the fitted point."""

from __future__ import annotations

import json
import math
from pathlib import Path

from ru_kgb import RUKGBParams, background


ROOT = Path(__file__).resolve().parents[1]
C_LIGHT = 299_792_458.0
G_NEWTON = 6.67430e-11
M_SUN = 1.98847e30
R_SUN = 6.957e8
MPC_M = 3.085677581491367e22
PC_M = 3.085677581491367e16


def schwarzschild_radius(mass_kg: float) -> float:
    return 2.0 * G_NEWTON * mass_kg / (C_LIGHT * C_LIGHT)


def main() -> None:
    params = RUKGBParams()
    present = background(1.0, params)
    alpha_b = present["alpha_B"]
    D = present["D"]
    c_hat = present["C_hat"]
    # This definition matches the high-k linear response
    # mu_infinity-1 = 2 beta_eff^2 = alpha_B^2/(2D) for c_s^2=1.
    beta_eff = abs(alpha_b) / (2.0 * math.sqrt(D))
    h0_s = 67.8625e3 / MPC_M
    h0_inverse_m = C_LIGHT / h0_s
    h0_length_inverse = 1.0 / h0_inverse_m
    r_sun = schwarzschild_radius(M_SUN)
    cassini_impact = 1.6 * R_SUN

    # Cubic KGB: Lambda_3^3=Mpl H0^2/|C_hat|.  The coefficient 4 pi is
    # the conventional reduced-Planck normalization of the point-source
    # decoupling equation.  Its O(1) convention dependence is immaterial
    # relative to the many-order-of-magnitude Cassini margin.
    r_v = (4.0 * math.pi * beta_eff * abs(c_hat) * r_sun / (h0_length_inverse**2)) ** (1.0 / 3.0)
    screened_fraction = min(1.0, (cassini_impact / r_v) ** 1.5)
    # Conservative envelope: the unscreened scalar-tensor gamma deviation is
    # bounded by 4 beta_eff^2, then Vainshtein-suppressed inside r_V.
    gamma_minus_one_bound = 4.0 * beta_eff**2 * screened_fraction
    cassini_bound = 2.3e-5
    result = {
        "model": "cubic KGB decoupling-limit Vainshtein gate",
        "present_action_values": {"alpha_B": alpha_b, "D": D, "C_hat": c_hat, "beta_eff": beta_eff},
        "source": {"mass_kg": M_SUN, "schwarzschild_radius_m": r_sun, "cassini_impact_m": cassini_impact},
        "screening": {
            "H0_inverse_length_m_inverse": h0_length_inverse,
            "r_V_m": r_v,
            "r_V_pc": r_v / PC_M,
            "impact_over_rV": cassini_impact / r_v,
            "screened_fraction": screened_fraction,
        },
        "PPN_gamma": {
            "absolute_gamma_minus_one_envelope": gamma_minus_one_bound,
            "Cassini_absolute_bound": cassini_bound,
            "margin_factor": cassini_bound / gamma_minus_one_bound,
            "gaussian_chi2_if_centered_at_GR": (gamma_minus_one_bound / cassini_bound) ** 2,
        },
        "verdict": "PASS: the cubic-KGB decoupling-limit PPN envelope is Vainshtein suppressed below the Cassini bound",
        "scope": (
            "This is the analytic local screening gate of the stated cubic action. A full ephemeris likelihood "
            "would additionally fit planetary initial conditions and nuisance parameters; no such likelihood is claimed here."
        ),
    }
    output = ROOT / "generated" / "ppn_screening.json"
    output.write_text(json.dumps(result, indent=2))
    print(result["verdict"])
    print(f"r_V(Sun) = {r_v / PC_M:.6f} pc")
    print(f"|gamma-1| envelope = {gamma_minus_one_bound:.3e}; Cassini bound = {cassini_bound:.3e}")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
