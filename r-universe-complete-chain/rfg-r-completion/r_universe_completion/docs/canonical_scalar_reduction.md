# Canonical-Scalar Quadratic Reduction of RFG-R

## Scope

This note completes one sharply defined piece of the scalar calculation: the
quadratic ADM reduction of RFG-R coupled to a single minimally coupled,
canonical scalar field.  We work in scalar-field comoving gauge,
`delta phi=0`, and set `H0=1` only in the executable code.  This calculation
does **not** replace the separately required dust, baryon, photon, neutrino,
and CDM perturbation systems.  It is a test of the sourced lapse--shift
constraints and their physical scalar mode.

## Intrinsic-curvature term

For `gamma_ij=a^2 exp(2 zeta) delta_ij` and `p=k^2/a^2`, the Fourier-space
curvature expansion is

```math
R3^(1)=4 p zeta,
R3^(2)=-10 p zeta^2.
```

Consequently, the terms in `N sqrt(gamma) Q R3` that contain neither lapse nor
`delta Q` include

```math
Q R3^(2)+3 Q zeta R3^(1)=2 Q p zeta^2.
```

This term is required.  Without it, the GR plus massless-canonical-scalar
limit produces `c_s^2=4/3` rather than the exact result `c_s^2=1`.

## Matter action and sourced constraints

For

```math
S_phi=int d^4x sqrt(-g)[-g^(mu nu) partial_mu phi partial_nu phi/2-U(phi)],
```

write `rho=P+U`, `P=phi_dot^2/2`, `pressure=P-U`, and
`m=phi_dot^2/M_Pl^2=(rho+pressure)/M_Pl^2`.  Dividing the quadratic density
by `(M_Pl^2/2)a^3`, the canonical field adds exactly

```math
Delta L_phi=m alpha^2-6(rho/M_Pl^2) alpha zeta
            +9(pressure/M_Pl^2) zeta^2.
```

It contains no `y=p beta` term in comoving gauge.  With `D=-F_XX/3`, the
combined nondynamical Hessian is

```math
H_(alpha,y) = [[2(-3 D H^2+m), 2 D H],
               [2 D H, 2(-D+2Q)/3]],

det H_(alpha,y)=-4 Delta/3,
Delta=6 D H^2 Q+(D-2Q)m.
```

Thus the constraints are invertible precisely when `Delta != 0`.

## Exact reduced action

Define `z=zeta_dot`, `x=(alpha,y)^T`, and

```math
h_z=(6 D H,-2D)^T,
h_zeta=(6 F-6F_X H-4 H Q_X p+4Qp-6 rho/M_Pl^2,
         4Q_X p/3)^T.
```

The complete quadratic density can be written as

```math
L=x^T H_(alpha,y) x/2+x^T(h_z z+h_zeta zeta)
  -3D z^2+(6F_X+4Q_X p) zeta z
  +9(pressure/M_Pl^2) zeta^2+2Qp zeta^2.
```

Eliminating `x` is an exact Schur complement.  It gives

```math
L_red=K zeta_dot^2+B(p) zeta zeta_dot+M(p)zeta^2,
K=-3D-h_z^T H_(alpha,y)^(-1) h_z/2
 =6 D Q m/[6 D H^2Q+(D-2Q)m].
```

Here `B=B0+B1 p` and `M=M0+M1 p+M2 p^2` are obtained without a
quasi-static approximation by the displayed matrix formula.  In particular,
the highest-gradient coefficient is

```math
M2=4[-3D Q^2+6H^2Q Q_X^2-12H Q^2 Q_X+6Q^3+Q_X^2 m]
   /(3 Delta),
```

in `H0=1` units.  Integrating the mixed term by parts gives

```math
C1=M1-H(B1_N+B1)/2,
C2=M2,
S_red=(M_Pl^2/2) int a^3 [K zeta_dot^2+(C0+C1p+C2p^2)zeta^2].
```

The physical low-`k` gradient coefficient is `G1=-C1`, the `k^4` coefficient
is `G2=-C2`, and the low-`k` sound speed is `c_s^2=G1/K` when `K>0`.

## Positivity on the reference branch

For the original high-`X` branch,

```math
D-2Q=theta(theta-1) Omega_R/(1+theta).
```

Therefore `theta=1.6`, `D>0`, `Q>0`, and canonical `m>0` imply
`Delta>0` and `K>0`.  The executable diagnostic uses a self-consistent,
flat reference background containing a massless canonical scalar with
`Omega_phi0=0.30` and `Omega_R0=0.70`.  It verifies on `0.03<=a<=1` that
`Delta`, `K`, `G1`, and `G2` are nonnegative (strictly positive except for
the GR-asymptotic `G2=0` limit).  The exact GR stiff-field limit is checked
analytically and in code: `K=6` and `c_s^2=1`.

Run it with:

```bash
python3 scripts/canonical_scalar.py
```

The output table is `generated/tables/canonical_scalar_reference.csv`.

## What This Does Not Establish

The calculation does not supply the multi-fluid kinetic and gradient matrix
needed for a physical matter/CMB universe, does not complete the full
metric--khronon Dirac count, and does not compute the cubic action or a
strong-coupling scale.  Those remain necessary before a Planck likelihood or
a comparison with LCDM can be claimed.
