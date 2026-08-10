# Exact photon--baryon--CDM--neutrino reduction

This document gives the action-level scalar reduction used by
`scripts/multifluid_reduction.py` and its independent audit in
`scripts/validate_multifluid_reduction.py`.  It removes the previously open
finite multi-fluid constraint calculation.  It does **not** report a CMB or
matter likelihood: no RFG-R chain has been run and no comparison with
LambdaCDM follows from this document.

The derivation has two layers which must not be conflated.

1. The lapse, shift, baryon, CDM, photon monopole/dipole and neutrino
   monopole/dipole variables have an exact finite quadratic action and exact
   algebraic Schur reduction.
2. Photon polarization and collisionless neutrino anisotropic stress are not
   perfect-fluid variables.  They are retained as the exact infinite
   Boltzmann hierarchy (or the unprojected massive-particle equation), rather
   than being terminated or replaced by a sound-speed ansatz.

The conventions are `H0=M_Pl=1`, `Delta -> -k^2`,
`s=k^2 psi/a^2`, and

```text
N=1+alpha,  N_i=partial_i psi,  gamma_ij=a^2 exp(2 zeta) delta_ij.
```

## 1. Starting action

The gravity action is the RFG-R ADM action mapped to the full extended EFT
basis, including `bar_m5 delta R^(3) delta K`.  For the scalar variables its
quadratic part is Eq. (85) of
[Frusciante, Papadomanolakis & Silvestri (2016)](https://arxiv.org/abs/1601.04064):

\[
\begin{aligned}
L_g={}&-W_0k^2\zeta^2-3a^2W_4\dot\zeta\alpha
-\frac32a^2W_5\dot\zeta^2-a^2W_4\alpha s-a^2W_5\dot\zeta s\\
&+a^4W_7s^2-2\bar m_5\frac{k^2}{a^2}\zeta s
+W_1\alpha^2-W_6k^2\alpha\zeta .
\end{aligned}
\]

Every `W_i(a)` is constructed from the analytic RFG-R action derivatives.
In particular, the implementation uses

\[
\dot{\bar m}_5=-\frac13 Q_{XX}\dot H,
\qquad
\bar m_5=-\frac13Q_X,
\]

and does not finite-difference either quantity when forming the constraint
matrix.

For a barotropic material component the exact Sorkin--Schutz action is

\[
S_i=-\int d^4x\,[\sqrt{-g}\,\rho_i(n_i)+J_i^\mu\partial_\mu\ell_i],
\qquad p_i=n_i\rho_{i,n}-\rho_i .
\]

Its scalar quadratic density, specialized only after `w_i=c_{s,i}^2` is
stated, is Eq. (III.12) of
[De Felice, Frusciante & Papadomanolakis (2017)](https://arxiv.org/abs/1609.03599):

\[
L_i=-\rho_i(1+w_i)(s+3\dot\zeta)v_i-\rho_i\dot\delta_i v_i
-\frac{\rho_i(1+w_i)}2\frac{k^2}{a^2}v_i^2
-\frac{\rho_iw_i}{2(1+w_i)}\delta_i^2-\rho_i\alpha\delta_i .
\]

For baryons and CDM, `w=0`.  For the photon and massless-neutrino **constraint
rows**, `w=1/3`.  This last statement does not close either kinetic species:
only its density and momentum moments enter the lapse/shift constraints.

## 2. Finite exact Schur reduction

Define

\[
\mathbf q=(\zeta,\delta_b,\delta_c,\delta_\gamma,\delta_\nu)^T,
\quad
\mathbf x=(\alpha,s,v_b,v_c,v_\gamma,v_\nu)^T.
\]

The complete finite density is

\[
L=\frac12\dot{\mathbf q}^{T}K_0\dot{\mathbf q}
+\mathbf x^TC\dot{\mathbf q}
+\frac12\mathbf x^TA\mathbf x
+\mathbf x^TD\mathbf q
+\frac12\mathbf q^TM_0\mathbf q.
\]

All nonzero matrix elements are:

\[
\begin{array}{lll}
(K_0)_{\zeta\zeta}=-3a^2W_5,&
C_{\alpha\zeta}=-3a^2W_4,&C_{s\zeta}=-a^2W_5,\\
A_{\alpha\alpha}=2W_1,&A_{\alpha s}=-a^2W_4,&A_{ss}=2a^4W_7,\\
D_{\alpha\zeta}=-W_6k^2,&D_{s\zeta}=-2\bar m_5k^2/a^2,&
(M_0)_{\zeta\zeta}=-2W_0k^2,\\
C_{v_i\zeta}=-3\rho_i(1+w_i),&C_{v_i\delta_i}=-\rho_i,&
A_{sv_i}=-\rho_i(1+w_i),\\
A_{v_iv_i}=-\rho_i(1+w_i)k^2/a^2,&D_{\alpha\delta_i}=-\rho_i,&
(M_0)_{\delta_i\delta_i}=-\rho_iw_i/(1+w_i).
\end{array}
\]

No term is inferred from a quasi-static relation, a `mu(a,k)` ansatz, or a
fluid closure.  Varying `alpha` and `s` gives the two sourced constraints

\[
\begin{aligned}
2W_1\alpha-a^2W_4s-3a^2W_4\dot\zeta-W_6k^2\zeta-\sum_i\rho_i\delta_i&=0,\\
-a^2W_4\alpha+2a^4W_7s-a^2W_5\dot\zeta
-2\bar m_5\frac{k^2}{a^2}\zeta-\sum_i\rho_i(1+w_i)v_i&=0.
\end{aligned}
\]

Their determinant is the exact identity

\[
\det A_{(\alpha,s)}=-a^4\left(W_4^2-4W_1W_7\right).
\]

Once all auxiliary variables are eliminated, the reduced matrices are

\[
K=K_0-C^TA^{-1}C,\qquad
B=-C^TA^{-1}D,\qquad
M=M_0-D^TA^{-1}D.
\]

This is the full finite algebraic reduction, including the nonzero `bar_m5`
operator.

## 3. Exact checks

`validate_multifluid_reduction.py` performs two distinct checks.

* **GR identity, exact arithmetic.**  At a rational GR background the code
  constructs the complete matrix with `fractions.Fraction`, takes the Schur
  complement exactly, and obtains `rank(K)=4` without a numerical tolerance.
  The finite monopole/dipole system therefore has four material scalar
  directions and no additional propagating gravity scalar.
* **RFG-R reference audit.**  The Planck 2018 species split below is supplied
  to the RFG-R reference branch on 25 by 17 `(a,k)` points.  The tested
  lapse--shift determinant residual is at most `4.148e-15`; after a diagonal
  congruence (which preserves inertia), every point has
  `(positive, negative, null)=(4,0,1)`.  The normalized null residual is at
  most `1.944e-13`.

The generated machine table is
[`multifluid_core_audit.csv`](../generated/tables/multifluid_core_audit.csv).
The numerical values test the analytic matrix identity; they are not fitted
RFG-R cosmological parameters and not a likelihood value.

## 4. Spatial metric equations and DAE closure

The spatially traceless scalar deformation must be retained until after
variation.  With

```text
gamma_ij=a^2 exp(2 zeta) exp(2 D_ij E),
D_ij=partial_i partial_j-delta_ij partial^2/3,
```

direct variation of the original ADM action and only then `E=0` gives, for
`s=k^2 beta/a^2`,

```math
dot(s) = -[3H+dot(Q)/Q+Q_X k^2/(3Qa^2)]s
         -(k^2/a^2)(alpha+zeta)
         +(Q_X/Q)(k^2/a^2)(H alpha-dot(zeta))
         -Pi/Q.
```

Here `Pi` is the total scalar anisotropic stress in `M_Pl^2 H0^2` units; a
massless species contributes `Pi=(4/5)rho Theta_2`.  The independently
derived GR identity is

```math
dot(s)=-3Hs-(k^2/a^2)(alpha+zeta)-Pi.
```

The complete scalar variation (lapse, shift, trace and traceless spatial
equations) is executable in
[`derive_spatial_traceless_equation.py`](../scripts/derive_spatial_traceless_equation.py).
It performs a direct Einstein-Hilbert ADM regression test for all four
residuals.

The finite Schur reduction has a precise null direction, not a missing
equation.  Define `Delta_i=delta_i+3(1+w_i)zeta`.  The transformed kinetic
matrix has a zero `zeta` row and column, and its symmetric mixing matrix gives
the algebraic curvature equation

```math
mu_zeta zeta + [M_(0i)-dot(B_(0i))] Delta_i=0,
mu_zeta=M_00-dot(B_00).
```

All rates are analytic action derivatives.  The closure audit finds that the
published RFG-R reference branch has `mu_zeta=0` at finite `(a,k)`; for
example `k/H0=2.51545672221` at `a=1`.  This is a linear DAE singularity, so
the reference branch does not define globally regular transfer functions.
The derivation, root table, and independent regression are in
[`rfg_dae_closure.md`](rfg_dae_closure.md),
[`rfg_dae_closure.py`](../scripts/rfg_dae_closure.py), and
[`multifluid_dae_closure.csv`](../generated/tables/multifluid_dae_closure.csv).

## 5. Photon and massless-neutrino hierarchy

The finite action is deliberately not used to eliminate `Theta_2` or any
higher moment.  In the same spatial gauge, the exact scalar intensity
hierarchy of [Hwang & Noh (2001)](https://arxiv.org/abs/astro-ph/0102005) is

\[
\dot\Theta_\ell=\frac{k}{a}\left[
\frac{\Theta_{\ell-1}}{2\ell-1}
-\frac{\Theta_{\ell+1}}{2\ell+3}\right]+M_\ell+C_\ell,
\]

with the `ell=0` streaming term equal to `-k Theta_1/(3a)` and

\[
M_0=-\dot\zeta-\frac{s}{3},\qquad
M_1=\frac{k}{a}\alpha,\qquad
M_2=\frac{2s}{3},\qquad M_{\ell\ge3}=0.
\]

The observable moment relations are

\[
\delta_\gamma=4\Theta_{\gamma0},\quad v_\gamma=\Theta_{\gamma1},\quad
\frac{\pi_\gamma}{\rho_\gamma}=\frac45\Theta_{\gamma2},
\]

and identically for a massless collisionless neutrino.  The `ell=0` equation
therefore reproduces the variational radiation continuity equation exactly;
the validation script checks this equality directly.

For photons, with `dot tau=n_e x_e sigma_T >= 0`,

\[
\begin{aligned}
C_0&=0,\\
C_1&=\dot\tau(v_b-\Theta_{\gamma1}),\\
C_2&=\dot\tau(P-\Theta_{\gamma2}),\\
C_{\ell\ge3}&=-\dot\tau\Theta_{\gamma\ell},\\
P&=\left(\Theta_{\gamma2}-\sqrt6E_2\right)/10.
\end{aligned}
\]

The scalar polarization hierarchy for `ell>=2` is

\[
\dot E_\ell=\frac{k}{a}\left[
\frac{\sqrt{\ell^2-4}}{2\ell-1}E_{\ell-1}
-\frac{\sqrt{(\ell+1)^2-4}}{2\ell+3}E_{\ell+1}\right]
-\dot\tau\left(E_\ell+\sqrt6P\,\delta_{\ell2}\right),
\]

with scalar `B_ell=0`.  For a collisionless massless neutrino, `C_ell=0` at
every `ell`.  These equations are defined for all `ell`; neither a tight
coupling approximation nor an `ell_max` closure enters this reduction.

The baryon/CDM dust equations in the same gauge are

\[
\dot\delta_{b,c}=-\frac{k}{a}v_{b,c}-3\dot\zeta-s,
\qquad
\dot v_c+Hv_c=\frac{k}{a}\alpha,
\]

and the baryon Euler equation receives the exactly momentum-conserving
Thomson term

\[
\dot v_b+Hv_b=\frac{k}{a}\alpha+
\frac{4\rho_\gamma}{3\rho_b}\dot\tau(v_\gamma-v_b)
\]

before the physical baryon pressure and recombination functions are supplied.
Those functions must come from a specified atomic-recombination calculation,
not a fitted effective sound speed.

## 6. Massive neutrinos

No massless approximation is made for a massive neutrino.  Its exact
collisionless phase-space equation uses
`epsilon(q,a)=sqrt(q^2+m_nu^2 a^2)`:

\[
\delta f'=-\frac{q}{\epsilon}\,\mathcal D_\gamma\delta f
+q\frac{\partial f_0}{\partial q}
\left[\frac{\epsilon}{q}A_{,\alpha}\gamma^\alpha
+(B_{\alpha|\beta}+C'_{\alpha\beta})\gamma^\alpha\gamma^\beta\right].
\]

The density, momentum and anisotropic stress entering the constraints and
evolution equations are the exact momentum integrals of this `delta f` with
the `epsilon` weights, as given in Eqs. (46)--(48) of Hwang & Noh.  The
function `massive_neutrino_delta_f_rhs` implements the local equation with
the full `epsilon`; it does not substitute `epsilon=q`.

## 7. Data provenance and scope

The 425-point matrix audit uses a published **reference record**, not an
RFG-R inference:

| quantity | value | provenance |
|---|---:|---|
| `h` | 0.6736 | Planck 2018 Table 2, base-LambdaCDM TT,TE,EE+lowE+lensing |
| `omega_b` | 0.02237 | same table |
| `omega_c` | 0.1200 | same table |
| `omega_gamma` | `2.4728e-5` | `T_CMB=2.7255 K` convention |
| `N_eff` | 3.046 | standard Planck baseline convention |

The code derives, rather than guesses,

\[
\omega_\nu=\omega_\gamma\frac78\left(\frac4{11}\right)^{4/3}N_{\rm eff},
\quad
\Omega_b=0.0493016923285,\quad
\Omega_c=0.264470410345,
\quad\Omega_r=9.21989275501\times10^{-5}.
\]

These are used only to exercise every matter row at a documented numerical
point.  They are conditional LambdaCDM posterior values and must never be
reported as a posterior, best fit, or likelihood result for RFG-R.

## Reproduction

```bash
cd r_universe_completion
PYTHONPATH=scripts python3 scripts/validate_multifluid_reduction.py
```

The complete package run is `bash scripts/run_all.sh`.
