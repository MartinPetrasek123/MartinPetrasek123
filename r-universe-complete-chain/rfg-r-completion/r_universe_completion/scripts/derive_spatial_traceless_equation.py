#!/usr/bin/env python3
"""Symbolically derive the scalar spatial RFG-R field equations.

The pre-existing scalar reduction fixed the spatial shear scalar to zero before
variation.  This script instead retains the traceless scalar deformation E,
expands the original ADM RFG-R action to second order for one Fourier mode,
and forms its Euler--Lagrange residuals.  The independent traceless (E) and
trace (zeta) spatial equations are both retained.  It is a derivation aid: no
spectrum or phenomenological closure is calculated here.

Conventions are N=1+alpha cos(kx), N_x=partial_x(beta cos(kx)), and
gamma_ij=a^2 exp[2 zeta cos(kx)] exp[2 D_ij(E cos(kx))], with
D_ij=partial_i partial_j-delta_ij partial^2/3.  The latter keeps E traceless.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class Symbols:
    """All amplitudes and background functions used by the symbolic expansion."""

    eps: sp.Symbol
    x: sp.Symbol
    a: sp.Symbol
    h: sp.Symbol
    h_dot: sp.Symbol
    h0: sp.Symbol
    k: sp.Symbol
    alpha: sp.Symbol
    alpha_dot: sp.Symbol
    zeta: sp.Symbol
    zeta_dot: sp.Symbol
    zeta_ddot: sp.Symbol
    shear: sp.Symbol
    shear_dot: sp.Symbol
    shear_ddot: sp.Symbol
    beta: sp.Symbol
    beta_dot: sp.Symbol
    q0: sp.Symbol
    q1: sp.Symbol
    q2: sp.Symbol
    q0_dot: sp.Symbol
    q1_dot: sp.Symbol
    q2_dot: sp.Symbol
    f0: sp.Symbol
    f1: sp.Symbol
    f2: sp.Symbol
    f0_dot: sp.Symbol
    f1_dot: sp.Symbol
    f2_dot: sp.Symbol
    eta: sp.Symbol
    sigma_completion: sp.Symbol
    geometry_completion: sp.Symbol


def symbols() -> Symbols:
    """Create real symbols; explicit names make the generated identity readable."""
    names = sp.symbols(
        "eps x a H Hdot H0 k alpha alphadot zeta zetadot zetaddot "
        "E Edot Eddot beta betadot Q QX QXX Qdot QXdot QXXdot "
        "F FX FXX Fdot FXdot FXXdot eta Sigma Xi",
        real=True,
    )
    return Symbols(*names)


def _order2(expression: sp.Expr, epsilon: sp.Symbol) -> sp.Expr:
    """Discard all terms cubic and higher in perturbation amplitudes."""
    return sp.expand(sp.series(sp.expand(expression), epsilon, 0, 3).removeO())


def _spatial_derivative(expression: sp.Expr, coordinate: sp.Symbol, epsilon: sp.Symbol) -> sp.Expr:
    return _order2(sp.diff(expression, coordinate), epsilon)


def _time_derivative(expression: sp.Expr, value: Symbols) -> sp.Expr:
    """Apply d/dt to a background/amplitude expression after mode averaging."""
    rates = {
        value.a: value.a * value.h,
        value.h: value.h_dot,
        value.alpha: value.alpha_dot,
        value.zeta: value.zeta_dot,
        value.zeta_dot: value.zeta_ddot,
        value.shear: value.shear_dot,
        value.shear_dot: value.shear_ddot,
        value.beta: value.beta_dot,
        value.q0: value.q0_dot,
        value.q1: value.q1_dot,
        value.q2: value.q2_dot,
        value.f0: value.f0_dot,
        value.f1: value.f1_dot,
        value.f2: value.f2_dot,
    }
    return sp.expand(sum(sp.diff(expression, variable) * rate for variable, rate in rates.items()))


def _mode_average(expression: sp.Expr, value: Symbols) -> sp.Expr:
    """Average a quadratic real Fourier mode over a spatial period exactly."""
    period = 2 * sp.pi / value.k
    averaged = sp.simplify(sp.integrate(sp.expand_trig(expression), (value.x, 0, period)) / period)
    # The physical scalar Fourier label is k>0. SymPy represents that branch
    # as the first entry of a Piecewise result when it integrates trig powers.
    if isinstance(averaged, sp.Piecewise):
        averaged = averaged.args[0][0]
    return sp.factor(averaged)


def quadratic_mode_lagrangian(value: Symbols | None = None) -> tuple[Symbols, sp.Expr]:
    """Return the mode-averaged O(eps^2) ADM gravitational density.

    The calculation retains the original trace--shear identity
    Q(X)[R3+sigma_ij sigma^ij]+2 H0^2 F(X), where X=K/(3H0).
    Q, QX, QXX and F, FX, FXX are independent background symbols at this
    stage; substituting the RFG-R identities is deliberately postponed.
    """
    value = value or symbols()
    e, x, a, k = value.eps, value.x, value.a, value.k
    cosine, sine = sp.cos(k * x), sp.sin(k * x)

    # D_xx E=-2 k^2 E/3 and D_yy E=D_zz E=k^2 E/3 for a mode along x.
    amplitude_x = value.zeta - 2 * k**2 * value.shear / 3
    amplitude_y = value.zeta + k**2 * value.shear / 3
    amplitude_x_dot = value.zeta_dot - 2 * k**2 * value.shear_dot / 3
    amplitude_y_dot = value.zeta_dot + k**2 * value.shear_dot / 3

    def metric_component(amplitude: sp.Expr) -> sp.Expr:
        return a**2 * (1 + 2 * e * amplitude * cosine + 2 * e**2 * amplitude**2 * cosine**2)

    gamma = [metric_component(amplitude_x), metric_component(amplitude_y), metric_component(amplitude_y)]
    # The inverse of exp(2 eps A cos(kx)) is known analytically to this order.
    # Writing it directly avoids a SymPy nseries failure on an expanded quotient.
    gamma_inverse = [
        a ** -2 * (1 - 2 * e * amplitude * cosine + 2 * e**2 * amplitude**2 * cosine**2)
        for amplitude in (amplitude_x, amplitude_y, amplitude_y)
    ]
    gamma_dot = [
        a**2
        * (
            2 * value.h * (1 + 2 * e * amplitude * cosine + 2 * e**2 * amplitude**2 * cosine**2)
            + 2 * e * amplitude_dot * cosine
            + 4 * e**2 * amplitude * amplitude_dot * cosine**2
        )
        for amplitude, amplitude_dot in ((amplitude_x, amplitude_x_dot), (amplitude_y, amplitude_y_dot), (amplitude_y, amplitude_y_dot))
    ]
    lapse = 1 + e * value.alpha * cosine
    lapse_inverse = _order2(1 / lapse, e)
    shift_covariant = [-e * k * value.beta * sine, sp.Integer(0), sp.Integer(0)]

    # Spatial Christoffel symbols for a diagonal metric that depends only on x.
    connection: list[list[list[sp.Expr]]] = [[[sp.Integer(0) for _ in range(3)] for _ in range(3)] for _ in range(3)]
    for upper in range(3):
        for lower_a in range(3):
            for lower_b in range(3):
                term = sp.Integer(0)
                if upper == lower_b:
                    term += _spatial_derivative(gamma[upper], x, e) if lower_a == 0 else 0
                if upper == lower_a:
                    term += _spatial_derivative(gamma[upper], x, e) if lower_b == 0 else 0
                if lower_a == lower_b and upper == 0:
                    term -= _spatial_derivative(gamma[lower_a], x, e)
                connection[upper][lower_a][lower_b] = _order2(gamma_inverse[upper] * term / 2, e)

    def covariant_shift_derivative(i: int, j: int) -> sp.Expr:
        partial = _spatial_derivative(shift_covariant[j], x, e) if i == 0 else sp.Integer(0)
        correction = sum(connection[upper][i][j] * shift_covariant[upper] for upper in range(3))
        return _order2(partial - correction, e)

    extrinsic_covariant = [
        _order2(lapse_inverse * (gamma_dot[i] - 2 * covariant_shift_derivative(i, i)) / 2, e)
        for i in range(3)
    ]
    extrinsic_mixed = [_order2(gamma_inverse[i] * extrinsic_covariant[i], e) for i in range(3)]
    trace_extrinsic = _order2(sum(extrinsic_mixed), e)
    shear_mixed = [_order2(component - trace_extrinsic / 3, e) for component in extrinsic_mixed]
    shear_square = _order2(sum(component * component for component in shear_mixed), e)

    # Build the three-dimensional Ricci scalar from its definition.
    ricci = [[sp.Integer(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            first = _spatial_derivative(connection[0][i][j], x, e)
            second = _spatial_derivative(sum(connection[upper][i][upper] for upper in range(3)), x, e) if j == 0 else 0
            third = sum(connection[upper][i][j] * connection[lower][upper][lower] for upper in range(3) for lower in range(3))
            fourth = sum(connection[lower][i][upper] * connection[upper][j][lower] for upper in range(3) for lower in range(3))
            ricci[i][j] = _order2(first - second + third - fourth, e)
    ricci_scalar = _order2(sum(gamma_inverse[index] * ricci[index][index] for index in range(3)), e)

    x_background = value.h / value.h0
    delta_x = _order2(trace_extrinsic / (3 * value.h0) - x_background, e)
    q = _order2(value.q0 + value.q1 * delta_x + value.q2 * delta_x**2 / 2, e)
    f = _order2(value.f0 + value.f1 * delta_x + value.f2 * delta_x**2 / 2, e)
    # sqrt(gamma) is independent of the traceless E amplitude exactly.  The
    # lapse is deliberately kept separate and multiplied only once below.
    measure = a**3 * (1 + 3 * e * value.zeta * cosine + sp.Rational(9, 2) * e**2 * value.zeta**2 * cosine**2)
    acceleration_square = _order2(
        gamma_inverse[0] * _spatial_derivative(sp.log(lapse), x, e) ** 2,
        e,
    )
    density = _order2(
        lapse
        * measure
        * (
            q * (ricci_scalar + shear_square)
            + value.sigma_completion * shear_square
            + value.geometry_completion * (ricci_scalar + shear_square)
            + 2 * value.h0**2 * f
            + value.eta * acceleration_square
        ),
        e,
    )
    quadratic = sp.expand(density).coeff(e, 2)
    return value, _mode_average(quadratic, value)


def spatial_traceless_residual(value: Symbols | None = None) -> tuple[Symbols, sp.Expr, sp.Expr]:
    """Return L^(2), dL/d(E_dot), and the exact Euler--Lagrange residual."""
    value, lagrangian = quadratic_mode_lagrangian(value)
    momentum = sp.factor(sp.diff(lagrangian, value.shear_dot))
    residual = sp.factor(_time_derivative(momentum, value) - sp.diff(lagrangian, value.shear))
    return value, sp.factor(lagrangian), residual


def spatial_trace_residual(value: Symbols | None = None) -> tuple[Symbols, sp.Expr, sp.Expr]:
    """Return L^(2), dL/d(zeta_dot), and the trace Euler--Lagrange residual.

    This variation is intentionally performed while the traceless scalar E is
    still present.  Fixing E=0 before the two spatial variations would erase
    the scalar anisotropic-stress equation and can conceal a constraint
    redundancy when one tries to evolve the metric.
    """
    value, lagrangian = quadratic_mode_lagrangian(value)
    momentum = sp.factor(sp.diff(lagrangian, value.zeta_dot))
    residual = sp.factor(_time_derivative(momentum, value) - sp.diff(lagrangian, value.zeta))
    return value, sp.factor(lagrangian), residual


def scalar_metric_residuals(value: Symbols | None = None) -> tuple[Symbols, dict[str, sp.Expr]]:
    """Return every scalar metric Euler residual before fixing E=0.

    ``alpha`` and ``beta`` are nondynamical and therefore give the Hamiltonian
    and momentum constraints.  ``zeta`` gives the spatial trace equation and
    ``E`` gives the spatial-traceless equation.  Keeping all four residuals
    together prevents an accidental replacement of the trace equation by a
    differentiated, dependent constraint.
    """
    value, lagrangian = quadratic_mode_lagrangian(value)
    zeta_momentum = sp.diff(lagrangian, value.zeta_dot)
    shear_momentum = sp.diff(lagrangian, value.shear_dot)
    return value, {
        "lapse": sp.factor(sp.diff(lagrangian, value.alpha)),
        "shift": sp.factor(sp.diff(lagrangian, value.beta)),
        "trace": sp.factor(_time_derivative(zeta_momentum, value) - sp.diff(lagrangian, value.zeta)),
        "traceless": sp.factor(_time_derivative(shear_momentum, value) - sp.diff(lagrangian, value.shear)),
    }


def spatial_gauge_metric_residuals(value: Symbols | None = None) -> tuple[Symbols, dict[str, sp.Expr]]:
    """Fix E=0 only after obtaining the four scalar metric residuals."""
    value, residuals = scalar_metric_residuals(value)
    fixed = {value.shear: 0, value.shear_dot: 0, value.shear_ddot: 0}
    return value, {name: sp.factor(residual.subs(fixed)) for name, residual in residuals.items()}


def spatial_gauge_shear_residual(value: Symbols | None = None) -> tuple[Symbols, sp.Expr]:
    """Return the E equation after varying first and then imposing E=0.

    The residual is normalized by the nonzero Fourier prefactor
    2 a k^4/3.  A matter anisotropic-stress variation is added only after the
    gravitational identity is independently checked below.
    """
    value, _, residual = spatial_traceless_residual(value)
    gauge_fixed = residual.subs({value.shear: 0, value.shear_dot: 0, value.shear_ddot: 0})
    normalized = sp.factor(gauge_fixed / (2 * value.a * value.k**4 / 3))
    return value, normalized


def spatial_gauge_trace_residual(value: Symbols | None = None) -> tuple[Symbols, sp.Expr]:
    """Return the trace spatial equation after the legitimate E=0 gauge fix."""
    value, _, residual = spatial_trace_residual(value)
    gauge_fixed = residual.subs({value.shear: 0, value.shear_dot: 0, value.shear_ddot: 0})
    return value, sp.factor(gauge_fixed)


def validate_gr_limit() -> None:
    """Verify all four variations against the Einstein--Hilbert ADM action."""
    value, residuals = spatial_gauge_metric_residuals()
    substitutions = {
        value.q0: 1,
        value.q1: 0,
        value.q2: 0,
        value.q0_dot: 0,
        value.q1_dot: 0,
        value.q2_dot: 0,
        # F=V-3X^2 for Einstein--Hilbert gravity plus a cosmological constant.
        value.f1: -6 * value.h / value.h0,
        value.f2: -6,
        value.f1_dot: -6 * value.h_dot / value.h0,
        value.f2_dot: 0,
        value.eta: 0,
        value.sigma_completion: 0,
        value.geometry_completion: 0,
    }
    gr = {name: sp.factor(residual.subs(substitutions)) for name, residual in residuals.items()}
    expected = {
        "lapse": value.a
        * (
            3 * value.f0 * value.h0**2 * value.a**2 * value.zeta
            - 6 * value.h**2 * value.a**2 * value.alpha
            + 18 * value.h**2 * value.a**2 * value.zeta
            + 6 * value.h * value.a**2 * value.zeta_dot
            + 2 * value.h * value.beta * value.k**2
            + 2 * value.k**2 * value.zeta
        ),
        "shift": 2 * value.a * value.k**2 * (value.h * value.alpha - value.zeta_dot),
        "trace": -value.a
        * (
            3 * value.f0 * value.h0**2 * value.a**2 * value.alpha
            + 9 * value.f0 * value.h0**2 * value.a**2 * value.zeta
            + 54 * value.h**2 * value.a**2 * value.zeta
            - 6 * value.h * value.a**2 * value.alpha_dot
            + 18 * value.h * value.a**2 * value.zeta_dot
            + 2 * value.h * value.beta * value.k**2
            - 6 * value.h_dot * value.a**2 * value.alpha
            + 18 * value.h_dot * value.a**2 * value.zeta
            + 6 * value.a**2 * value.zeta_ddot
            + 2 * value.alpha * value.k**2
            + 2 * value.beta_dot * value.k**2
            + 2 * value.k**2 * value.zeta
        ),
        "traceless": 2 * value.a * value.k**4 * (-value.beta_dot - value.h * value.beta - value.alpha - value.zeta) / 3,
    }
    for name, identity in expected.items():
        if sp.simplify(gr[name] - identity) != 0:
            raise AssertionError(f"incorrect GR {name} limit: {sp.sstr(gr[name])}")


def main() -> None:
    value, lagrangian, residual = spatial_traceless_residual()
    _, gauge_fixed = spatial_gauge_shear_residual(value)
    _, _, trace_residual = spatial_trace_residual(value)
    _, trace_gauge_fixed = spatial_gauge_trace_residual(value)
    _, metric_gauge_fixed = spatial_gauge_metric_residuals(value)
    validate_gr_limit()
    print("Mode-averaged L^(2):")
    print(sp.sstr(lagrangian))
    print("\nSpatial-traceless E Euler--Lagrange residual:")
    print(sp.sstr(residual))
    print("\nLaTex residual:")
    print(sp.latex(residual))
    print("\nE=0 gauge shear residual (GR limit checked):")
    print(sp.sstr(gauge_fixed))
    print("\nSpatial-trace zeta Euler--Lagrange residual:")
    print(sp.sstr(trace_residual))
    print("\nE=0 gauge trace residual:")
    print(sp.sstr(trace_gauge_fixed))
    print("\nE=0 gauge lapse and shift residuals:")
    print(f"lapse = {sp.sstr(metric_gauge_fixed['lapse'])}")
    print(f"shift = {sp.sstr(metric_gauge_fixed['shift'])}")


if __name__ == "__main__":
    main()
