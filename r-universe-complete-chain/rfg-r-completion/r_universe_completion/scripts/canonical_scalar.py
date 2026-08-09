#!/usr/bin/env python3
"""Canonical-scalar quadratic reduction for the regular RFG-R action.

This module is deliberately limited to one minimally coupled, canonical scalar
field in comoving gauge (delta phi = 0).  It is not a replacement for the
dust, radiation, neutrino, or photon-baryon perturbation systems required by
a Boltzmann calculation.

Units: H0=1.  The reduced action is normalized as

  S^(2) = (M_Pl^2/2) int dt d^3k a^3 [K zeta_dot^2 + C(k) zeta^2].

The intrinsic-curvature contribution +2 Q p zeta^2 is included explicitly.
It is essential for the General-Relativity limit.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import math
from pathlib import Path

from rfg_regularized import (
    RFGRegularizedParams,
    Q,
    Q_prime,
    potential,
    response,
    response_prime,
)


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class MasslessScalarReference:
    """Flat RFG-R branch sourced by a canonical massless scalar.

    The source is Omega_phi0 a^-6.  Its canonical scalar has rho=P and
    m=phi_dot^2/M_Pl^2=2 rho/M_Pl^2.  ``omega_m0`` below fixes only the
    present response amplitude Omega_R0=1-Omega_phi0; no dust is present in
    this diagnostic background.
    """

    omega_phi0: float = 0.30
    theta: float = 1.6
    epsilon: float = 1.0e-8
    p: int = 4

    @property
    def parameters(self) -> RFGRegularizedParams:
        return RFGRegularizedParams(
            omega_m0=self.omega_phi0,
            omega_r0=0.0,
            theta=self.theta,
            epsilon=self.epsilon,
            p=self.p,
        )

    def source(self, a: float) -> float:
        return self.omega_phi0 * a**-6

    def solve_E(self, a: float) -> float:
        """Solve E^2=Omega_phi0 a^-6+R_epsilon(E) on the positive branch."""
        params = self.parameters
        lo = 0.0
        hi = max(1.0, math.sqrt(self.source(a)) + 1.0)

        def residual(e: float) -> float:
            return e * e - self.source(a) - response(e, params)

        while residual(hi) <= 0.0:
            hi *= 2.0
        for _ in range(320):
            mid = 0.5 * (lo + hi)
            if residual(mid) > 0.0:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)


def _second_derivative(function, x: float, params: RFGRegularizedParams) -> float:
    step = max(1.0e-6 * x, 1.0e-7)
    left = max(0.0, x - step)
    return (function(x + step, params) - function(left, params)) / (x + step - left)


def _dot_inverse(u: tuple[float, float], inverse, v: tuple[float, float]) -> float:
    return u[0] * (inverse[0][0] * v[0] + inverse[0][1] * v[1]) + u[1] * (
        inverse[1][0] * v[0] + inverse[1][1] * v[1]
    )


def _kernel_polynomials(
    h: float,
    d: float,
    q: float,
    q_x: float,
    f: float,
    f_x: float,
    rho: float,
    pressure: float,
    m: float,
) -> dict[str, float]:
    """Eliminate lapse and shift by an exact two-by-two Schur complement."""
    a_alpha = 2.0 * (-3.0 * d * h * h + m)
    a_alpha_y = 2.0 * d * h
    a_y = 2.0 * (-d + 2.0 * q) / 3.0
    determinant = a_alpha * a_y - a_alpha_y * a_alpha_y
    inverse = (
        (a_y / determinant, -a_alpha_y / determinant),
        (-a_alpha_y / determinant, a_alpha / determinant),
    )
    h_z = (6.0 * d * h, -2.0 * d)
    kinetic = -3.0 * d - 0.5 * _dot_inverse(h_z, inverse, h_z)

    def at_p(wave_p: float) -> tuple[float, float]:
        h_zeta = (
            6.0 * f - 6.0 * f_x * h - 4.0 * h * q_x * wave_p + 4.0 * q * wave_p - 6.0 * rho,
            4.0 * q_x * wave_p / 3.0,
        )
        mix = 6.0 * f_x + 4.0 * q_x * wave_p - _dot_inverse(h_z, inverse, h_zeta)
        mass = 9.0 * pressure - 0.5 * _dot_inverse(h_zeta, inverse, h_zeta)
        # From Q R^(2) + 3 Q zeta R^(1), see the companion derivation.
        mass += 2.0 * q * wave_p
        return mix, mass

    b0, m0 = at_p(0.0)
    b_at_1, m_at_1 = at_p(1.0)
    b_at_2, m_at_2 = at_p(2.0)
    b1 = b_at_1 - b0
    m2 = 0.5 * (m_at_2 - 2.0 * m_at_1 + m0)
    m1 = m_at_1 - m0 - m2
    return {
        "K": kinetic,
        "B0": b0,
        "B1": b1,
        "M0": m0,
        "M1": m1,
        "M2": m2,
        "constraint_determinant": determinant,
    }


def instantaneous_reduction(a: float, background: MasslessScalarReference) -> dict[str, float]:
    """Return all quadratic coefficients at one scale factor.

    The function evaluates the exact RFG-R functions, reconstructs
    F=-3X^2Q+V and D=-F_XX/3, and then eliminates (alpha,y).
    """
    params = background.parameters
    h = background.solve_E(a)
    q = Q(h, params)
    q_x = Q_prime(h, params)
    q_xx = _second_derivative(Q_prime, h, params)
    source_f = h * h - response(h, params) - h * h * q - h**3 * q_x
    source_f_x = (
        2.0 * h
        - response_prime(h, params)
        - 2.0 * h * q
        - 4.0 * h * h * q_x
        - h**3 * q_xx
    )
    v = potential(h, params)
    v_x = (v - 3.0 * source_f) / h
    v_xx = -3.0 * source_f_x / h
    f = -3.0 * h * h * q + v
    f_x = -6.0 * h * q - 3.0 * h * h * q_x + v_x
    f_xx = -6.0 * q - 12.0 * h * q_x - 3.0 * h * h * q_xx + v_xx
    d = -f_xx / 3.0

    # rho/M_Pl^2=3 Omega_phi0 a^-6 in H0=1 units for a canonical massless field.
    rho = 3.0 * background.source(a)
    pressure = rho
    m = 2.0 * rho
    values = _kernel_polynomials(h, d, q, q_x, f, f_x, rho, pressure, m)
    delta = 6.0 * d * h * h * q + (d - 2.0 * q) * m
    values.update(
        {
            "a": a,
            "E": h,
            "Q": q,
            "Q_X": q_x,
            "D": d,
            "rho_over_Mpl2": rho,
            "pressure_over_Mpl2": pressure,
            "m": m,
            "Delta": delta,
            "background_residual": h * h - background.source(a) - response(h, params),
            "background_relative_residual": (h * h - background.source(a) - response(h, params)) / (h * h),
        }
    )
    return values


def reduced_coefficients(a: float, background: MasslessScalarReference) -> dict[str, float]:
    """Add the integration-by-parts gradient coefficients C1 and C2.

    With B=B0+B1 p and p=k^2/a^2, integration by parts gives
    C1=M1-H(B1_N+B1)/2 and C2=M2.  The action has
    C=C0+C1 p+C2 p^2, hence G1=-C1 and G2=-C2.
    """
    values = instantaneous_reduction(a, background)
    step = 2.0e-5
    plus = instantaneous_reduction(a * math.exp(step), background)
    minus = instantaneous_reduction(a * math.exp(-step), background)
    b1_n = (plus["B1"] - minus["B1"]) / (2.0 * step)
    c1 = values["M1"] - 0.5 * values["E"] * (b1_n + values["B1"])
    values.update(
        {
            "B1_N": b1_n,
            "C1": c1,
            "C2": values["M2"],
            "G1": -c1,
            "G2": -values["M2"],
            "c_s2_low_k": -c1 / values["K"],
        }
    )
    return values


def validate_gr_limit() -> None:
    """Exact GR + massless-scalar check of the curvature term and normalization."""
    # On the GR stiff-field solution, D=2, Q=1, Q_X=0, m=6H^2 and rho=P=3H^2.
    # The Schur complement gives K=6, B1=4/H and M1=2.  Since H_N=-3H,
    # C1=2-H[(12/H)+(4/H)]/2=-6, i.e. c_s^2=(-C1)/K=1.
    h = 7.0
    values = _kernel_polynomials(
        h=h,
        d=2.0,
        q=1.0,
        q_x=0.0,
        f=-3.0 * h * h,
        f_x=-6.0 * h,
        rho=3.0 * h * h,
        pressure=3.0 * h * h,
        m=6.0 * h * h,
    )
    assert abs(values["K"] - 6.0) < 1.0e-12
    assert abs(values["B1"] - 4.0 / h) < 1.0e-12
    assert abs(values["M1"] - 2.0) < 1.0e-12
    c1 = values["M1"] - 0.5 * h * (12.0 / h + values["B1"])
    assert abs(c1 + 6.0) < 1.0e-12


def main() -> None:
    validate_gr_limit()
    background = MasslessScalarReference()
    # At a<0.03 the RFG-R correction is below numerical double-precision
    # resolution; the exact GR-limit test above covers that asymptotic regime.
    scales = [0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
    rows = [reduced_coefficients(a, background) for a in scales]
    for row in rows:
        assert abs(row["background_relative_residual"]) < 2.0e-12
        assert row["Delta"] > 0.0
        assert row["K"] > 0.0
        assert row["G1"] > 0.0
        assert row["G2"] >= -1.0e-9

    output = ROOT / "generated" / "tables" / "canonical_scalar_reference.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "a", "E", "D", "Q", "Delta", "K", "G1", "G2", "c_s2_low_k", "constraint_determinant"
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows([{field: row[field] for field in fields} for row in rows])

    present = rows[-1]
    print("Canonical-scalar RFG-R reduction validation OK")
    print("exact GR stiff-field limit: K=6, c_s^2=1")
    print(f"min Delta on reference grid = {min(row['Delta'] for row in rows):.10g}")
    print(f"min K on reference grid     = {min(row['K'] for row in rows):.10g}")
    print(f"min G1 on reference grid    = {min(row['G1'] for row in rows):.10g}")
    print(f"min G2 on reference grid    = {min(row['G2'] for row in rows):.10g}")
    print(f"a=1 low-k scalar c_s^2     = {present['c_s2_low_k']:.10g}")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
