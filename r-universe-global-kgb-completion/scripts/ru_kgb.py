#!/usr/bin/env python3
"""Globally completed, luminal KGB realization of the fitted R-alpha branch.

The physical background equals the R-alpha equation used in the late-time
likelihood for every observed redshift (0 < a <= 1):

  E^2 = Omega_m a^-3 + Omega_r a^-4 + Omega_R E^(alpha a).

For a > 1 the exponent is smoothly saturated.  This does not change a single
past-light-cone datum, but it makes the covariant action single-valued for all
phi and gives a future de Sitter limit.  The action is

  S = integral sqrt(-g) [Mpl^2 R/2 + A(phi) X + B(phi) X^2
                         - V(phi) - C(phi) X Box(phi)] + S_m,

where phi/Mpl = ln(a) on the reconstructed FLRW solution.  A, B, C and V are
defined functions of phi through the unique positive background root, rather
than an interpolation table that only exists on one numerical trajectory.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class RUKGBParams:
    omega_m0: float = 0.3022734375
    omega_r0: float = 9.083909e-5
    alpha: float = 0.497421875
    a_saturation: float = 2.0
    target_cs2: float = 1.0

    @property
    def omega_R0(self) -> float:
        return 1.0 - self.omega_m0 - self.omega_r0

    @property
    def beta_radiation(self) -> float:
        """Fixed, non-fitted radiation regularizer for the scalar kinetic term.

        Its magnitude is determined entirely by the already fitted R-alpha
        coefficient and the radiation density; it is not an additional sampled
        parameter.  The physical early-time alpha_B floor uses the opposite
        sign, which is required for a positive scalar gradient numerator in
        radiation domination.  For alpha -> 0 it vanishes and the exact LCDM
        limit is recovered.
        """
        return self.alpha * self.omega_r0

    def validate(self) -> None:
        if not (0.0 <= self.omega_m0 < 1.0):
            raise ValueError("omega_m0 must lie in [0, 1)")
        if not (0.0 <= self.omega_r0 < 1.0):
            raise ValueError("omega_r0 must lie in [0, 1)")
        if not (0.0 < self.omega_R0 < 1.0):
            raise ValueError("omega_R0 must be positive")
        if not (0.0 <= self.alpha < 1.0):
            raise ValueError("alpha must lie in [0, 1) on the stable branch")
        if not (self.a_saturation > 1.0):
            raise ValueError("a_saturation must exceed one")
        if not (self.alpha * self.a_saturation < 2.0):
            raise ValueError("the saturated response exponent must remain below two")
        if not (self.target_cs2 > 0.0):
            raise ValueError("target_cs2 must be positive")


def _smoothstep(t: float) -> tuple[float, float]:
    """C-infinity step and its derivative on 0 < t < 1."""
    if t <= 0.0:
        return 0.0, 0.0
    if t >= 1.0:
        return 1.0, 0.0
    left = math.exp(-1.0 / t)
    right = math.exp(-1.0 / (1.0 - t))
    denom = left + right
    value = left / denom
    # d ln(left/right)/dt = 1/t^2 + 1/(1-t)^2.
    derivative = value * (1.0 - value) * (1.0 / (t * t) + 1.0 / ((1.0 - t) ** 2))
    return value, derivative


def a_effective(a: float, params: RUKGBParams) -> tuple[float, float]:
    """Return a_eff and d a_eff/dN for the global finite-window completion."""
    params.validate()
    if a <= 0.0:
        raise ValueError("a must be positive")
    a_sat = params.a_saturation
    if a <= 1.0:
        return a, a
    if a >= a_sat:
        return a_sat, 0.0
    t = (a - 1.0) / (a_sat - 1.0)
    step, step_a = _smoothstep(t)
    step_a /= a_sat - 1.0
    effective = a + (a_sat - a) * step
    effective_a = 1.0 - step + (a_sat - a) * step_a
    return effective, a * effective_a


def matter_sources(a: float, params: RUKGBParams) -> tuple[float, float]:
    return params.omega_m0 * a**-3, params.omega_r0 * a**-4


def _response(e: float, exponent: float, params: RUKGBParams) -> float:
    return params.omega_R0 * e**exponent


def solve_E(a: float, params: RUKGBParams, iterations: int = 260) -> float:
    """Unique expanding root of the globally completed R-alpha equation."""
    params.validate()
    effective, _ = a_effective(a, params)
    exponent = params.alpha * effective
    matter, radiation = matter_sources(a, params)

    def residual(e: float) -> float:
        return e * e - _response(e, exponent, params) - matter - radiation

    lo = 0.0
    hi = max(1.0, math.sqrt(matter + radiation + params.omega_R0) + 1.0)
    while residual(hi) <= 0.0:
        hi *= 2.0
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def background(a: float, params: RUKGBParams) -> dict[str, float]:
    """Evaluate the exact KGB reconstruction in Mpl=H0=1 units."""
    params.validate()
    e = solve_E(a, params)
    effective, effective_n = a_effective(a, params)
    exponent = params.alpha * effective
    exponent_n = params.alpha * effective_n
    matter, radiation = matter_sources(a, params)
    r = _response(e, exponent, params)
    f_e = 2.0 * e - r * exponent / e
    e_n = (r * exponent_n * math.log(e) - 3.0 * matter - 4.0 * radiation) / f_e
    r_n = r * (exponent_n * math.log(e) + exponent * e_n / e)

    omega_m = matter / (e * e)
    omega_r = radiation / (e * e)
    omega_R = r / (e * e)

    # The original R-alpha braiding target is retained for a <= 1.  A fixed
    # negative radiation floor removes the vanishing-Qs branch and makes the
    # standard Horndeski scalar gradient numerator positive in radiation
    # domination.  Its magnitude adds no cosmological parameter.
    raw_gate = exponent / (1.0 + exponent)
    raw_gate_n = exponent_n / ((1.0 + exponent) ** 2)
    omega_R_n = omega_R * (r_n / r - 2.0 * e_n / e)
    beta_floor = params.beta_radiation
    alpha_b = -beta_floor + (1.0 + beta_floor) * raw_gate * omega_R
    alpha_b_n = (1.0 + beta_floor) * (raw_gate_n * omega_R + raw_gate * omega_R_n)

    # Bellini-Sawicki scalar sound-speed numerator for alpha_M=alpha_T=0.
    # c_s^2=N_s/D and D=alpha_K+3 alpha_B^2/2.  The matter/radiation
    # enthalpy is the final term.  This is the convention associated with
    # Q_s=2 Mpl^2 D/(2-alpha_B)^2.
    sound_numerator = (
        (2.0 - alpha_b) * (-e_n / e + 0.5 * alpha_b)
        + alpha_b_n
        - (3.0 * omega_m + 4.0 * omega_r)
    )
    if sound_numerator <= 0.0:
        raise RuntimeError("chosen braiding branch has a scalar gradient instability")

    D = sound_numerator / params.target_cs2
    alpha_k = D - 1.5 * alpha_b * alpha_b
    if D <= 0.0:
        raise RuntimeError("scalar kinetic coefficient is not positive")
    Q_s = 2.0 * D / ((2.0 - alpha_b) ** 2)
    F_s = Q_s * sound_numerator / D
    c_s2 = sound_numerator / D

    # The functions below define the covariant action off shell as functions
    # of phi/Mpl=N, because a=exp(phi/Mpl) and E(phi) is the positive root.
    X = 0.5 * e * e
    C = alpha_b / (e * e)
    C_n = (alpha_b_n - 2.0 * alpha_b * e_n / e) / (e * e)
    rho_phi = 3.0 * r
    # This is algebraically identical to the acceleration equation expression
    # -(2 E E_N+3 E^2)-Omega_r a^-4, but is stable in the radiation era.
    p_phi = -3.0 * r - r_n
    R_rho = rho_phi - 6.0 * e * e * X * C + 2.0 * C_n * X * X
    R_p = p_phi + 2.0 * C_n * X * X + 2.0 * C * X * (e * e_n)
    B = (alpha_k * e * e - (R_rho + R_p) + 8.0 * C_n * X * X - 12.0 * e * e * X * C) / (8.0 * X * X)
    A = (R_rho + R_p - 4.0 * B * X * X) / (2.0 * X)
    V = 0.5 * (R_rho - R_p) - B * X * X

    rho_reconstructed = A * X + 3.0 * B * X * X + V + 6.0 * e * e * X * C - 2.0 * C_n * X * X
    p_reconstructed = A * X + B * X * X - V - 2.0 * C_n * X * X - 2.0 * C * X * (e * e_n)
    background_residual = e * e - matter - radiation - r
    # rho_phi,N + 3(rho_phi+p_phi) is the homogeneous scalar equation
    # multiplied by phi_dot.  Keep this direct form; rearranging it through
    # the Friedmann equation needlessly subtracts radiation-era quantities.
    conservation_residual = 3.0 * r_n + 3.0 * (rho_phi + p_phi)

    return {
        "a": a,
        "z": 1.0 / a - 1.0,
        "phi_over_Mpl": math.log(a),
        "a_effective": effective,
        "E": e,
        "E_N": e_n,
        "response_exponent": exponent,
        "Omega_m": omega_m,
        "Omega_r": omega_r,
        "Omega_R": omega_R,
        "rho_R": r,
        "rho_R_N": r_n,
        "w_R": -1.0 - r_n / (3.0 * r),
        "alpha_B": alpha_b,
        "alpha_B_N": alpha_b_n,
        "alpha_K": alpha_k,
        "D": D,
        "sound_numerator": sound_numerator,
        "F_s": F_s,
        "Q_s_over_Mpl2": Q_s,
        "c_s2": c_s2,
        "A_hat": A,
        "B_hat": B,
        "C_hat": C,
        "C_hat_phi": C_n,
        "V_hat": V,
        "rho_phi_target": rho_phi,
        "p_phi_target": p_phi,
        "rho_phi_reconstructed": rho_reconstructed,
        "p_phi_reconstructed": p_reconstructed,
        "Q_T_over_Mpl2": 1.0,
        "c_T2": 1.0,
        "dL_GW_over_dL_EM": 1.0,
        "background_residual": background_residual,
        "conservation_residual": conservation_residual,
        "rho_reconstruction_residual": rho_reconstructed - rho_phi,
        "p_reconstruction_residual": p_reconstructed - p_phi,
    }


def trajectory(
    params: RUKGBParams,
    a_min: float = 1.0e-8,
    a_max: float = 1.0e3,
    points: int = 3001,
) -> list[dict[str, float]]:
    if not (0.0 < a_min < a_max):
        raise ValueError("require 0 < a_min < a_max")
    ns = np.linspace(math.log(a_min), math.log(a_max), points)
    return [background(float(math.exp(n)), params) for n in ns]
