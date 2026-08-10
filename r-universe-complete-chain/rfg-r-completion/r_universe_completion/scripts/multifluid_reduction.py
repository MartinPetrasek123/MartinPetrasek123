#!/usr/bin/env python3
"""Exact photon--baryon--CDM--neutrino reduction for the RFG-R action.

The finite quadratic action contains only the monopole/dipole variables needed
by the lapse and shift constraints.  Photons and collisionless neutrinos are
*not* closed as perfect fluids: their exact scalar Boltzmann hierarchy is
given below for every ell >= 0.  The massive-neutrino equation is retained in
its unprojected phase-space form, so no massless-neutrino or fluid closure is
hidden in this module.

Conventions: H0=M_Pl=1, spatially flat FLRW, Fourier Delta -> -k^2, and
s=k^2 psi/a^2.  The ADM perturbations are
N=1+alpha, N_i=partial_i psi, gamma_ij=a^2 exp(2 zeta) delta_ij.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np

from extended_eft_mapping import extended_eft_coefficients
from extended_eft_scalar_stability import _w_coefficients
from rfg_regularized import RFGRegularizedParams


@dataclass(frozen=True)
class FluidSpecies:
    """A barotropic monopole/dipole source used only in the constraint action."""

    name: str
    omega0: float
    w: float

    def rho(self, a: float) -> float:
        """Return rho/M_Pl^2 H0^2 on the stated homogeneous background."""
        return 3.0 * self.omega0 * a ** (-3.0 * (1.0 + self.w))


@dataclass(frozen=True)
class Planck2018Reference:
    """Published base-LambdaCDM reference point, never an RFG-R fit.

    Values are the Planck 2018 TT,TE,EE+lowE+lensing posterior means in
    Table 2 of Aghanim et al., A&A 641, A6 (2020), with the standard fixed
    T_CMB and N_eff convention used in that table.  It is used only to give
    a data-recorded species split for the internal matrix audit.
    """

    h: float = 0.6736
    omega_b_h2: float = 0.02237
    omega_c_h2: float = 0.1200
    omega_gamma_h2: float = 2.4728e-5
    n_eff: float = 3.046

    @property
    def omega_nu_h2(self) -> float:
        return self.omega_gamma_h2 * (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * self.n_eff

    def species(self) -> tuple[FluidSpecies, ...]:
        h2 = self.h * self.h
        return (
            FluidSpecies("baryon", self.omega_b_h2 / h2, 0.0),
            FluidSpecies("cdm", self.omega_c_h2 / h2, 0.0),
            FluidSpecies("photon_monopole_dipole", self.omega_gamma_h2 / h2, 1.0 / 3.0),
            FluidSpecies("massless_neutrino_monopole_dipole", self.omega_nu_h2 / h2, 1.0 / 3.0),
        )

    def rfg_background_parameters(self, *, theta: float = 1.6) -> RFGRegularizedParams:
        """Use the data point only as an audit input, not as an RFG-R posterior."""
        h2 = self.h * self.h
        return RFGRegularizedParams(
            omega_m0=(self.omega_b_h2 + self.omega_c_h2) / h2,
            omega_r0=(self.omega_gamma_h2 + self.omega_nu_h2) / h2,
            theta=theta,
        )


@dataclass(frozen=True)
class ScalarMetricSource:
    """Metric variables entering the kinetic equations in the RFG spatial gauge."""

    alpha: float
    zeta_dot: float
    s: float


def quadratic_blocks(
    a: float,
    k_over_h0: float,
    params: RFGRegularizedParams,
    species: Iterable[FluidSpecies],
) -> dict[str, object]:
    """Return the exact finite monopole/dipole quadratic action blocks.

    With q=(zeta,delta_i) and x=(alpha,s,v_i), the scalar density is

      L = qdot^T K0 qdot/2 + x^T C qdot + x^T A x/2
          + x^T D q + q^T M0 q/2.

    It follows by adding the Sorkin--Schutz action for each barotropic source
    to the full extended-EFT gravity action (including bar_m5 deltaR3 deltaK).
    For photon and neutrino rows it supplies exactly the ell=0,1 constraint
    sector only; their ell>=2 hierarchy remains separate below.
    """
    if a <= 0.0 or k_over_h0 <= 0.0:
        raise ValueError("a and k_over_h0 must be positive")
    sources = tuple(species)
    values = _w_coefficients(a, params)
    w0, w1 = values["W0"], values["W1"]
    w4, w5, w6, w7 = values["W4"], values["W5"], values["W6"], values["W7"]
    m5 = values["m5_bar_hat"]
    n = 1 + len(sources)
    m = 2 + len(sources)
    k2 = k_over_h0 * k_over_h0
    p = k2 / (a * a)
    K0 = np.zeros((n, n), dtype=float)
    C = np.zeros((m, n), dtype=float)
    A = np.zeros((m, m), dtype=float)
    D = np.zeros((m, n), dtype=float)
    M0 = np.zeros((n, n), dtype=float)

    # Gravity: Frusciante, Papadomanolakis & Silvestri, arXiv:1601.04064,
    # Eqs. (85)--(86), expressed through s=k^2 psi/a^2.
    K0[0, 0] = -3.0 * a * a * w5
    C[0, 0] = -3.0 * a * a * w4
    C[1, 0] = -a * a * w5
    A[0, 0] = 2.0 * w1
    A[0, 1] = A[1, 0] = -a * a * w4
    A[1, 1] = 2.0 * a**4 * w7
    D[0, 0] = -w6 * k2
    D[1, 0] = -2.0 * m5 * p
    M0[0, 0] = -2.0 * w0 * k2

    # Matter: De Felice, Frusciante & Papadomanolakis, arXiv:1609.03599,
    # Eq. (III.12), specialized only after the equation of state is stated.
    for number, source in enumerate(sources):
        q_index = 1 + number
        x_index = 2 + number
        rho = source.rho(a)
        one_plus_w = 1.0 + source.w
        C[x_index, 0] = -3.0 * rho * one_plus_w
        C[x_index, q_index] = -rho
        A[1, x_index] = A[x_index, 1] = -rho * one_plus_w
        A[x_index, x_index] = -rho * one_plus_w * p
        D[0, q_index] = -rho
        M0[q_index, q_index] = -rho * source.w / one_plus_w

    return {
        "K0": K0,
        "C": C,
        "A": A,
        "D": D,
        "M0": M0,
        "a": a,
        "k_over_H0": k_over_h0,
        "p": p,
        "species": sources,
        "constraint_discriminant": w4 * w4 - 4.0 * w1 * w7,
        "constraint_block_determinant": float(np.linalg.det(A[:2, :2])),
        **values,
    }


def reduce_auxiliaries(blocks: dict[str, object]) -> dict[str, np.ndarray]:
    """Eliminate all nondynamical x variables by an exact Schur complement."""
    K0 = np.asarray(blocks["K0"], dtype=float)
    C = np.asarray(blocks["C"], dtype=float)
    A = np.asarray(blocks["A"], dtype=float)
    D = np.asarray(blocks["D"], dtype=float)
    M0 = np.asarray(blocks["M0"], dtype=float)
    solved_c = np.linalg.solve(A, C)
    solved_d = np.linalg.solve(A, D)
    return {
        "K": K0 - C.T @ solved_c,
        "B": -C.T @ solved_d,
        "M": M0 - D.T @ solved_d,
    }


def solve_lapse_shift(
    blocks: dict[str, object],
    zeta_dot: float,
    deltas_dot: Iterable[float],
    zeta: float,
    deltas: Iterable[float],
    velocities: Iterable[float],
) -> tuple[float, float]:
    """Solve the sourced alpha,s constraint pair without eliminating v_i."""
    C = np.asarray(blocks["C"], dtype=float)
    A = np.asarray(blocks["A"], dtype=float)
    D = np.asarray(blocks["D"], dtype=float)
    qdot = np.array((zeta_dot, *deltas_dot), dtype=float)
    q = np.array((zeta, *deltas), dtype=float)
    v = np.array(tuple(velocities), dtype=float)
    if len(v) != A.shape[0] - 2:
        raise ValueError("one velocity is required for every fluid source")
    rhs = -(C[:2] @ qdot + D[:2] @ q + A[:2, 2:] @ v)
    alpha, s = np.linalg.solve(A[:2, :2], rhs)
    return float(alpha), float(s)


def temperature_metric_source(ell: int, metric: ScalarMetricSource, k_over_a: float) -> float:
    """Exact scalar metric source M_ell in the RFG spatial gauge.

    Hwang--Noh use chi=-psi in this convention, hence M_0=-zeta_dot-s/3,
    M_1=(k/a)alpha, and M_2=2s/3 for K=0.
    """
    if ell == 0:
        return -metric.zeta_dot - metric.s / 3.0
    if ell == 1:
        return k_over_a * metric.alpha
    if ell == 2:
        return 2.0 * metric.s / 3.0
    return 0.0


def massless_hierarchy_rhs(
    ell: int,
    theta_lower: float,
    theta: float,
    theta_upper: float,
    metric: ScalarMetricSource,
    k_over_a: float,
    collision: float = 0.0,
) -> float:
    """One exact ell equation for a photon or massless-neutrino hierarchy.

    ``collision`` is the already projected C_ell.  Calling this local equation
    for every ell>=0 defines an infinite hierarchy; no l_max or terminal
    closure is part of this function.
    """
    if ell < 0:
        raise ValueError("ell must be non-negative")
    if ell == 0:
        streaming = -k_over_a * theta_upper / 3.0
    else:
        streaming = k_over_a * (theta_lower / (2.0 * ell - 1.0) - theta_upper / (2.0 * ell + 3.0))
    return streaming + temperature_metric_source(ell, metric, k_over_a) + collision


def photon_collision(
    ell: int,
    theta_ell: float,
    theta_2: float,
    e_2: float,
    baryon_velocity: float,
    optical_depth_dot: float,
) -> float:
    """Exact Thomson collision multipole for scalar photon perturbations."""
    if optical_depth_dot < 0.0:
        raise ValueError("this convention uses optical_depth_dot=n_e x_e sigma_T >= 0")
    polarization_source = (theta_2 - math.sqrt(6.0) * e_2) / 10.0
    if ell == 0:
        return 0.0
    if ell == 1:
        return optical_depth_dot * (baryon_velocity - theta_ell)
    if ell == 2:
        return optical_depth_dot * (polarization_source - theta_ell)
    return -optical_depth_dot * theta_ell


def photon_e_polarization_rhs(
    ell: int,
    e_lower: float,
    e_ell: float,
    e_upper: float,
    theta_2: float,
    e_2: float,
    k_over_a: float,
    optical_depth_dot: float,
) -> float:
    """Exact scalar E-polarization hierarchy (ell>=2; scalar B_ell=0)."""
    if ell < 2:
        raise ValueError("scalar E polarization starts at ell=2")
    kappa_lower = math.sqrt(float(ell * ell - 4))
    kappa_upper = math.sqrt(float((ell + 1) * (ell + 1) - 4))
    streaming = k_over_a * (kappa_lower * e_lower / (2.0 * ell - 1.0) - kappa_upper * e_upper / (2.0 * ell + 3.0))
    polarization_source = (theta_2 - math.sqrt(6.0) * e_2) / 10.0
    return streaming - optical_depth_dot * (e_ell + math.sqrt(6.0) * polarization_source * (ell == 2))


def baryon_cdm_continuity_rhs(zeta_dot: float, s: float, velocity: float, k_over_a: float) -> float:
    """Exact dust continuity equation in the same gauge as the action."""
    return -k_over_a * velocity - 3.0 * zeta_dot - s


def dust_euler_rhs(alpha: float, velocity: float, hubble: float, k_over_a: float, momentum_exchange: float = 0.0) -> float:
    """Exact dust Euler equation; momentum_exchange is J/(rho) in cosmic time."""
    return -hubble * velocity + k_over_a * alpha + momentum_exchange


def photon_baryon_momentum_exchange(
    rho_gamma: float, rho_baryon: float, optical_depth_dot: float, photon_velocity: float, baryon_velocity: float
) -> float:
    """The Thomson term in the baryon Euler equation from total momentum conservation."""
    return 4.0 * rho_gamma * optical_depth_dot * (photon_velocity - baryon_velocity) / (3.0 * rho_baryon)


def massless_anisotropic_stress(rho: float, theta_2: float) -> float:
    """Return Pi=(4/5) rho Theta_2 in the scalar kinetic convention."""
    return 4.0 * rho * theta_2 / 5.0


def spatial_traceless_shear_rhs_from_coefficients(
    *,
    a: float,
    k_over_h0: float,
    hubble: float,
    q: float,
    q_x: float,
    q_dot: float,
    alpha: float,
    zeta: float,
    zeta_dot: float,
    shear: float,
    anisotropic_stress: float,
) -> float:
    """Return dot(s) from the scalar spatial-traceless field equation.

    This is obtained by retaining the traceless scalar spatial metric mode E,
    varying the original ADM RFG-R action, and only then fixing E=0.  The
    derivation is reproduced in ``derive_spatial_traceless_equation.py``.
    Units are H0=M_Pl=1 and ``s=k^2 beta/a^2``. ``anisotropic_stress`` is the
    total Pi/M_Pl^2 H0^2 of the photon, massless-neutrino and massive-neutrino
    kinetic moments.
    """
    if a <= 0.0 or k_over_h0 <= 0.0 or q == 0.0:
        raise ValueError("a, k_over_h0 and q must be nonzero on the physical branch")
    gradient = k_over_h0 * k_over_h0 / (a * a)
    return (
        -(3.0 * hubble + q_dot / q + q_x * gradient / (3.0 * q)) * shear
        - gradient * (alpha + zeta)
        + gradient * hubble * q_x * alpha / q
        - gradient * q_x * zeta_dot / q
        - anisotropic_stress / q
    )


def spatial_traceless_shear_rhs(
    a: float,
    k_over_h0: float,
    params: RFGRegularizedParams,
    *,
    alpha: float,
    zeta: float,
    zeta_dot: float,
    shear: float,
    anisotropic_stress: float,
) -> float:
    """Evaluate the action-derived RFG-R shear equation on the background."""
    row = extended_eft_coefficients(a, params)
    h = float(row["E"])
    q = float(row["Q"])
    q_x = float(row["Q_X"])
    q_dot = q_x * float(row["Hdot_over_H0_sq"])
    return spatial_traceless_shear_rhs_from_coefficients(
        a=a,
        k_over_h0=k_over_h0,
        hubble=h,
        q=q,
        q_x=q_x,
        q_dot=q_dot,
        alpha=alpha,
        zeta=zeta,
        zeta_dot=zeta_dot,
        shear=shear,
        anisotropic_stress=anisotropic_stress,
    )


def massive_neutrino_delta_f_rhs(
    q: float,
    a: float,
    mass: float,
    directional_streaming: float,
    alpha_gradient_projection: float,
    shear_projection: float,
    df0_dq: float,
) -> float:
    """Linear collisionless massive-neutrino equation for delta f in conformal time.

    This is Eq. (43) of Hwang & Noh, arXiv:astro-ph/0102005, after separating
    f=f0(q)+delta f.  The caller supplies the Fourier/angular projections and
    df0_dq.  epsilon=sqrt(q^2+m^2 a^2) is retained exactly; the equation has
    no fluid or massless approximation.
    """
    if q <= 0.0 or a <= 0.0 or mass < 0.0:
        raise ValueError("q,a must be positive and mass non-negative")
    epsilon = math.sqrt(q * q + mass * mass * a * a)
    redshift = (epsilon / q) * alpha_gradient_projection + shear_projection
    return -(q / epsilon) * directional_streaming + redshift * q * df0_dq


def normalized_inertia(matrix: np.ndarray, *, rtol: float = 2.0e-10) -> tuple[int, int, int, float]:
    """Return inertia after a diagonal congruence, preserving all exact signs."""
    diagonal = np.abs(np.diag(matrix))
    if np.any(diagonal <= 0.0):
        raise ValueError("kinetic diagonal must be nonzero before normalization")
    scale = 1.0 / np.sqrt(diagonal)
    normalized = scale[:, None] * matrix * scale[None, :]
    eigenvalues = np.linalg.eigvalsh(normalized)
    threshold = rtol * max(1.0, float(np.max(np.abs(eigenvalues))))
    positive = int(np.count_nonzero(eigenvalues > threshold))
    negative = int(np.count_nonzero(eigenvalues < -threshold))
    null = len(eigenvalues) - positive - negative
    null_residual = float(np.min(np.abs(eigenvalues)) / max(1.0, np.max(np.abs(eigenvalues))))
    return positive, negative, null, null_residual
