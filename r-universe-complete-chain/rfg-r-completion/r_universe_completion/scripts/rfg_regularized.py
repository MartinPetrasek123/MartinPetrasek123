#!/usr/bin/env python3
"""Regularized R-Universe / RFG-R background and EFT coefficient functions.

The completion has two jobs:
1. retain the reconstructed RFG branch for X much larger than epsilon;
2. recover the Einstein-Hilbert ADM action at X=0.

All quantities are dimensionless apart from the explicit factors of H0 in the
paper.  This module intentionally has no external numerical dependency.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class RFGRegularizedParams:
    omega_m0: float = 0.30000
    omega_r0: float = 9.0e-5
    theta: float = 1.6
    # Published conditional-branch representative.  The separate
    # regularization audit demonstrates why this is not a globally healthy
    # completion and must not be treated as one.
    epsilon: float = 1.0e-8
    p: int = 4

    @property
    def omega_R0(self) -> float:
        return 1.0 - self.omega_m0 - self.omega_r0

    @property
    def A(self) -> float:
        return self.omega_R0 / (1.0 + self.theta)

    def validate(self) -> None:
        if not (0.0 < self.theta < 2.0):
            raise ValueError("theta must be in (0, 2)")
        if self.p <= self.theta or self.p % 2 != 0:
            raise ValueError("p must be an even integer larger than theta")
        if not (0.0 < self.epsilon < 1.0):
            raise ValueError("epsilon must be in (0, 1)")
        if not (0.0 < self.omega_R0 < 1.0):
            raise ValueError("Omega_R0 must be in (0, 1)")


def _positive_x(x: float) -> float:
    if x < 0.0:
        raise ValueError("the cosmological branch uses X >= 0")
    return x


def _denominator(x: float, p: RFGRegularizedParams) -> float:
    return x**p.p + p.epsilon**p.p


def _transition_terms(x: float, p: RFGRegularizedParams) -> tuple[float, float]:
    """Return s=X^p/(X^p+epsilon^p) and A*f without overflow.

    The regular response can be evaluated at very early cosmological times only
    after factoring out the larger of X and epsilon.  In particular, the
    direct X**p representation overflows for the high-order regularizations
    needed to keep Q positive throughout the transition.
    """
    if x == 0.0:
        return 0.0, 0.0
    nu = 1.0 + p.theta / p.p
    if x >= p.epsilon:
        ratio = (p.epsilon / x) ** p.p
        denom = 1.0 + ratio
        s = 1.0 / denom
        correction = p.A * x ** (-p.theta) / denom**nu
    else:
        ratio = (x / p.epsilon) ** p.p
        denom = 1.0 + ratio
        s = ratio / denom
        correction = p.A * p.epsilon ** (-p.theta) * ratio / denom**nu
    return s, correction


def response(x: float, p: RFGRegularizedParams) -> float:
    """Regularized relational term in E^2 = S + response(E)."""
    _positive_x(x)
    if x == 0.0:
        return 0.0
    nu = 1.0 + p.theta / p.p
    if x >= p.epsilon:
        ratio = (p.epsilon / x) ** p.p
        return p.omega_R0 * x ** (2.0 - p.theta) / (1.0 + ratio) ** nu
    ratio = (x / p.epsilon) ** p.p
    return (
        p.omega_R0
        * p.epsilon ** (2.0 - p.theta)
        * (x / p.epsilon) ** (p.p + 2)
        / (1.0 + ratio) ** nu
    )


def response_prime(x: float, p: RFGRegularizedParams) -> float:
    """d response / dX, evaluated analytically."""
    _positive_x(x)
    if x == 0.0:
        return 0.0
    nu = 1.0 + p.theta / p.p
    s, _ = _transition_terms(x, p)
    return response(x, p) * (p.p + 2.0 - nu * p.p * s) / x


def Q(x: float, p: RFGRegularizedParams) -> float:
    """Regularized coefficient multiplying intrinsic and extrinsic curvature."""
    _positive_x(x)
    _, correction = _transition_terms(x, p)
    return 1.0 - correction


def Q_prime(x: float, p: RFGRegularizedParams) -> float:
    _positive_x(x)
    if x == 0.0:
        return 0.0
    nu = 1.0 + p.theta / p.p
    s, correction = _transition_terms(x, p)
    return -correction * p.p * (1.0 - nu * s) / x


def Q_second(x: float, p: RFGRegularizedParams) -> float:
    """Stable central difference for the coefficient table and checks."""
    if x == 0.0:
        return 0.0 if p.p > 2 else math.nan
    step = max(1.0e-7 * x, 1.0e-14)
    return (Q_prime(x + step, p) - Q_prime(max(0.0, x - step), p)) / (x + step - max(0.0, x - step))


def Q_minimum(p: RFGRegularizedParams) -> tuple[float, float]:
    """Return the exact global minimum of Q and the X at which it occurs."""
    ratio = p.p / p.theta
    x_star = p.epsilon * ratio ** (1.0 / p.p)
    correction_max = (
        p.A
        * p.epsilon ** (-p.theta)
        * ratio
        / (1.0 + ratio) ** (1.0 + p.theta / p.p)
    )
    return 1.0 - correction_max, x_star


def original_response(x: float, p: RFGRegularizedParams) -> float:
    _positive_x(x)
    return p.omega_R0 * x ** (2.0 - p.theta)


def original_Q(x: float, p: RFGRegularizedParams) -> float:
    _positive_x(x)
    return 1.0 - p.A * x ** (-p.theta)


def source_for_potential(x: float, p: RFGRegularizedParams) -> float:
    """F in V - X V_X = 3 F for the exact regularized branch."""
    # Evaluate X^2(1-Q) directly.  Writing X^2-X^2Q loses the leading
    # low-X term to cancellation when Q rounds to 1 in double precision.
    _, correction = _transition_terms(x, p)
    return x * x * correction - response(x, p) - x**3 * Q_prime(x, p)


def potential(x: float, p: RFGRegularizedParams, points: int | None = None) -> float:
    """The representative with no linear-in-X boundary contribution.

    V(X) = -3 X integral_0^X F(s)/s^2 ds.
    The regularization makes the integrand O(s^p) at the lower endpoint.
    """
    _positive_x(x)
    if x == 0.0:
        return 0.0
    if points is None:
        # Resolve increasingly sharp but smooth high-p transitions.  The
        # p=4 legacy representative retains 4097 points; p=32 uses 32769.
        refinement = max(1, math.ceil(p.p / 4.0))
        points = 4097 * refinement - (refinement - 1)
    # The transition is centered at epsilon, so a linear grid would miss it by
    # many orders of magnitude.  Integrate in log X after an analytic origin cap.
    origin_cap = p.epsilon * 1.0e-4
    if x <= origin_cap:
        return potential_small_x_coefficient(p) * x ** (p.p + 2)
    f_small = p.omega_R0 * (p.p - p.theta) / ((1.0 + p.theta) * p.epsilon ** (p.p + p.theta))
    integral_origin = f_small * origin_cap ** (p.p + 1) / (p.p + 1.0)
    log_grid = np.linspace(math.log(origin_cap), math.log(x), points)
    grid = np.exp(log_grid)
    nu = 1.0 + p.theta / p.p
    r = np.empty_like(grid)
    q = np.empty_like(grid)
    qx = np.empty_like(grid)
    q_correction = np.empty_like(grid)
    upper = grid >= p.epsilon
    lower = ~upper
    if np.any(upper):
        ratio = (p.epsilon / grid[upper]) ** p.p
        denom = 1.0 + ratio
        s = 1.0 / denom
        r[upper] = p.omega_R0 * grid[upper] ** (2.0 - p.theta) / denom**nu
        correction = p.A * grid[upper] ** (-p.theta) / denom**nu
        q[upper] = 1.0 - correction
        q_correction[upper] = correction
        qx[upper] = -correction * p.p * (1.0 - nu * s) / grid[upper]
    if np.any(lower):
        ratio = (grid[lower] / p.epsilon) ** p.p
        denom = 1.0 + ratio
        s = ratio / denom
        r[lower] = (
            p.omega_R0
            * p.epsilon ** (2.0 - p.theta)
            * (grid[lower] / p.epsilon) ** (p.p + 2)
            / denom**nu
        )
        correction = p.A * p.epsilon ** (-p.theta) * ratio / denom**nu
        q[lower] = 1.0 - correction
        q_correction[lower] = correction
        qx[lower] = -correction * p.p * (1.0 - nu * s) / grid[lower]
    f_values = grid * grid * q_correction - r - grid**3 * qx
    integral = integral_origin + float(np.trapz((f_values / (grid * grid)) * grid, log_grid))
    return -3.0 * x * integral


def potential_prime(x: float, p: RFGRegularizedParams) -> float:
    if x == 0.0:
        return 0.0
    return (potential(x, p) - 3.0 * source_for_potential(x, p)) / x


def potential_small_x_coefficient(p: RFGRegularizedParams) -> float:
    """Coefficient C in V = C X^(p+2) + higher orders near X=0."""
    return -3.0 * p.omega_R0 * (p.p - p.theta) / (
        (1.0 + p.theta) * (p.p + 1.0) * p.epsilon ** (p.p + p.theta)
    )


def background_source(a: float, p: RFGRegularizedParams) -> float:
    return p.omega_m0 * a**-3 + p.omega_r0 * a**-4


def background_residual(E: float, a: float, p: RFGRegularizedParams) -> float:
    return E * E - background_source(a, p) - response(E, p)


def solve_E(a: float, p: RFGRegularizedParams, tol: float = 1.0e-13) -> float:
    """Positive branch root by bracketed bisection."""
    p.validate()
    if not a > 0.0:
        raise ValueError("a must be positive")
    lo = 0.0
    hi = max(1.0, math.sqrt(background_source(a, p)) + 1.0)
    while background_residual(hi, a, p) <= 0.0:
        hi *= 2.0
    for _ in range(320):
        mid = 0.5 * (lo + hi)
        value = background_residual(mid, a, p)
        if abs(value) < tol:
            return mid
        if value > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def observables(a: float, p: RFGRegularizedParams) -> dict[str, float]:
    E = solve_E(a, p)
    source = background_source(a, p)
    matter = p.omega_m0 * a**-3 / (E * E)
    radiation = p.omega_r0 * a**-4 / (E * E)
    relational = response(E, p) / (E * E)
    dE_dN = (-3.0 * p.omega_m0 * a**-3 - 4.0 * p.omega_r0 * a**-4) / (2.0 * E - response_prime(E, p))
    dlnE_dN = dE_dN / E
    q_t = Q(E, p)
    alpha_m = Q_prime(E, p) * dE_dN / q_t
    return {
        "a": a,
        "z": 1.0 / a - 1.0,
        "E": E,
        "Omega_m": matter,
        "Omega_r": radiation,
        "Omega_R": relational,
        "closure": matter + radiation + relational,
        "dlnE_dln_a": dlnE_dN,
        "q": -1.0 - dlnE_dN,
        "w_R": -1.0 - (response_prime(E, p) * dE_dN / response(E, p)) / 3.0,
        "Q_T": q_t,
        "alpha_M": alpha_m,
        "dL_GW_over_dL_EM": math.sqrt(Q(1.0, p) / q_t),
        "source": source,
    }


def find_z_acc(p: RFGRegularizedParams) -> float:
    lo, hi = 0.0, 10.0
    q_lo = observables(1.0 / (1.0 + lo), p)["q"]
    q_hi = observables(1.0 / (1.0 + hi), p)["q"]
    if q_lo * q_hi > 0.0:
        raise RuntimeError("acceleration transition not bracketed")
    for _ in range(150):
        mid = 0.5 * (lo + hi)
        q_mid = observables(1.0 / (1.0 + mid), p)["q"]
        if q_lo * q_mid <= 0.0:
            hi = mid
        else:
            lo = mid
            q_lo = q_mid
    return 0.5 * (lo + hi)


def eft_coefficients(a: float, p: RFGRegularizedParams) -> dict[str, float]:
    """Independent ADM derivatives needed by a full 3+1 EFT implementation.

    The entries use X=E(a), R^(3)=0 and K_ij=H gamma_ij on FLRW.  Dimensionful
    H0 factors are restored by the implementation from the formulas in the paper.
    """
    row = observables(a, p)
    x = row["E"]
    q = Q(x, p)
    qx = Q_prime(x, p)
    qxx = Q_second(x, p)
    v = potential(x, p)
    vx = potential_prime(x, p)
    # V_XX by a stable finite difference. The term is used only as a table value;
    # the defining equation for V_X is exact above.
    step = max(1.0e-6 * x, 1.0e-10)
    vxx = (potential_prime(x + step, p) - potential_prime(max(0.0, x - step), p)) / (
        x + step - max(0.0, x - step)
    )
    row.update(
        {
            "X": x,
            "Q_X": qx,
            "Q_XX": qxx,
            "V": v,
            "V_X": vx,
            "V_XX": vxx,
            "L_R3": q,
            "L_KR3_trace": qx / 3.0,
            "F_XX": -6.0 * q - 12.0 * x * qx - 3.0 * x * x * qxx + vxx,
        }
    )
    return row
