# Exact RFG-R Multi-Fluid DAE Closure And Singular Constraint Surface

This note closes the finite scalar action reduction algebraically. It does not
insert a quasi-static relation, choose a temporal gauge, truncate a
photon/neutrino hierarchy, or finite-difference an action coefficient.

The result is decisive for the published regular RFG-R reference branch: the
scalar curvature constraint is singular on a continuous surface in `(a,k)`.
Consequently that branch cannot supply globally regular linear CMB, lensing,
or matter-transfer functions, and no RFG-R likelihood may be reported for it.

## 1. Finite action and the null direction

The exact finite action from
[`photon_baryon_cdm_neutrino_reduction.md`](photon_baryon_cdm_neutrino_reduction.md)
is

```math
L = dot(q)^T K_0 dot(q)/2 + x^T C dot(q) + x^T A x/2
  + x^T D q + q^T M_0 q/2,
```

with `q=(zeta,delta_b,delta_c,delta_gamma,delta_nu)^T` and the lapse, shift,
and velocity auxiliaries in `x`. Define exactly

```math
Delta_i = delta_i + 3(1+w_i) zeta,
q = T y,
y=(zeta,Delta_b,Delta_c,Delta_gamma,Delta_nu)^T,
T_(i0)=-3(1+w_i), T_(ii)=1.
```

After the exact Schur complement,

```math
K = K_0-C^T A^(-1)C,
B = -C^T A^(-1)D,
M = M_0-D^T A^(-1)D,
```

the transformed density is

```math
L_red = dot(y)^T (T^TKT) dot(y)/2 + dot(y)^T(T^TBT)y
      + y^T(T^TMT)y/2.
```

The action-derived identities are

```math
(T^T K T)_(0j)=0,                 B=B^T.
```

Thus `zeta` has no kinetic term. This is not a time-gauge choice: preferred
foliation time is already fixed by the action. It is an algebraic constraint
of the multi-fluid action.

## 2. Exact curvature constraint

Integrating the `dot(zeta)` term by parts gives

```math
mu_zeta(a,k) zeta + c_i(a,k) Delta_i
 + u_i(a,k) dot(Delta_i)=0,
```

where

```math
mu_zeta = M_00-dot(B_00),
c_i = M_(0i)-dot(B_(0i)),
u_i = B_(i0)-B_(0i).
```

For the RFG-R action, `u_i=0` identically because the reduced `B` matrix is
symmetric. Away from a zero of `mu_zeta`, the action itself fixes

```math
zeta = -c_i Delta_i / mu_zeta.
```

The coefficients and their time derivatives are calculated analytically. In
particular, the implementation evaluates `Q_X`, `Q_XX`, `Q_XXX`, `V_XX`,
`V_XXX`, `W_i`, and `dot(W_i)` from the regularized action, and differentiates
the Schur complement as

```math
dot(B) = -dot(C)^T A^(-1)D
         + C^T A^(-1)dot(A)A^(-1)D
         - C^T A^(-1)dot(D).
```

No derivative of a tabulated background is used in the calculation.

## 3. Reproducible result

Run:

```bash
cd r_universe_completion
PYTHONPATH=scripts python3 scripts/rfg_dae_closure.py
PYTHONPATH=scripts python3 scripts/validate_rfg_dae_closure.py
```

The first command writes
[`multifluid_dae_closure.csv`](../generated/tables/multifluid_dae_closure.csv).
On 49 logarithmic scale factors from `a=1e-7` to `a=1`, the calculation gives

```text
max transformed zeta kinetic residual = 5.923e-15
max reduced B antisymmetry residual    = 1.360e-16
curvature-constraint roots on grid     = 24 / 49
a=1 root k/H0                          = 2.51545672221
```

The independent validator checks the analytic `dot(B)` expression against a
separate centered numerical derivative only as a regression test; its maximum
relative residual is `1.114e-10`. It also brackets the `a=1` zero with a
normalized residual below `1e-8`.

Examples from the generated table are

| a | z | root k/H0 |
|---:|---:|---:|
| 0.0348070 | 27.7298 | 325.888461 |
| 0.0953162 | 9.4914 | 58.7660194 |
| 0.261016 | 2.8312 | 11.0993570 |
| 0.510897 | 0.9573 | 4.26905912 |
| 1 | 0 | 2.51545672 |

At the root, the equation that normally determines `zeta` loses rank. This is
a singularity of the exact linear differential-algebraic system, not an
incomplete Boltzmann implementation and not a numerical hierarchy cutoff.

## 4. Consequence

Photon and neutrino equations remain exact infinite kinetic hierarchies, and
the spatial-traceless metric equation including their anisotropic stress is
derived in [`derive_spatial_traceless_equation.py`](../scripts/derive_spatial_traceless_equation.py).
Those facts cannot remove a zero of `mu_zeta`: before a transfer function can
cross that surface, the finite metric constraint itself ceases to determine
the curvature perturbation.

Therefore the reference RFG-R branch is excluded as a globally regular
linear-theory candidate on the stated domain. A different action would have
to be specified and re-derived before a CMB/matter/PPN likelihood could be
meaningfully evaluated. The executed GR CAMB and Planck calculations remain
GR interface references only and are not evidence for RFG-R.
