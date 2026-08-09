# RFG-R: Complete Definition And Derivations

## 1. Purpose And Scope

The original RFG branch fixes a background response proportional to
`X^(2-theta)` and a tensor normalization `Q=1-A X^-theta`. For `theta>0`, the
latter is singular at `X=0`, which prevents a standard asymptotically flat
weak-field expansion from being a statement of the action. RFG-R resolves that
specific defect while retaining the branch where it was derived: `X=H/H0` of
order unity or larger.

RFG-R is a low-energy multiscale EFT. It has a cosmological RFG-R action and a
locally matched GR domain. The matching rule is part of the model definition;
it is not an unstated screening assumption.

## 2. Regularized Cosmological Action

Let `p` be an even integer with `p>theta`; the reference choice is `p=4`.
Define

```math
nu=1+theta/p,
```

```math
R_epsilon(X)=Omega_R0 X^(p+2)/(X^p+epsilon^p)^nu,
```

```math
Q_epsilon(X)=1-A X^p/(X^p+epsilon^p)^nu,
```

```math
A=Omega_R0/(1+theta).
```

The action is

```math
S_cos=(M_Pl^2/2) integral d^4x sqrt(-g)
 [Q_epsilon(X)(R3+K_mn K^mn-K^2)+2H0^2 V_epsilon(X)]+S_m[g,Psi].
```

Here `X=K/(3H0)`, the khronon determines the hypersurfaces, and all matter
fields `Psi` couple minimally and universally to `g_mn`.

## 3. Exact Background Reconstruction

The desired regularized FLRW equation is

```math
X^2=rho/(3 M_Pl^2 H0^2)+R_epsilon(X).
```

The lapse variation of the action gives

```math
X^2 Q_epsilon+X^3 Q_epsilon,X+(V_epsilon-XV_epsilon,X)/3
=rho/(3 M_Pl^2 H0^2).
```

Therefore the potential is fixed by

```math
V_epsilon-XV_epsilon,X=3F_epsilon(X),
```

```math
F_epsilon(X)=X^2-R_epsilon(X)-X^2Q_epsilon(X)-X^3Q_epsilon,X(X).
```

Writing `V_epsilon=X Y_epsilon` gives `Y_epsilon,X=-3F_epsilon/X^2` and hence

```math
V_epsilon(X)=-3X integral_0^X ds F_epsilon(s)/s^2.
```

The integrand is finite at the origin. Any extra `C X` is the familiar ADM
boundary representative and changes none of the field equations.

## 4. Recovery Of The Original Branch

For `X >> epsilon`,

```math
R_epsilon(X)=Omega_R0 X^(2-theta)[1-(1+theta/p)(epsilon/X)^p+O((epsilon/X)^(2p))],
```

```math
Q_epsilon(X)=1-A X^-theta[1-(1+theta/p)(epsilon/X)^p+O((epsilon/X)^(2p))].
```

Thus `epsilon=1e-8` and `p=4` give an analytic relative deformation below
`2.5e-32` for `X>=0.8`, which covers the original expanding cosmological
branch. Double precision verification is necessarily limited to about `1e-15`.

## 5. Local Einstein-Hilbert Limit

At small `X`,

```math
Q_epsilon(X)=1-A epsilon^(-p-theta)X^p+O(X^(2p)),
```

```math
V_epsilon(X)=
-[3 Omega_R0(p-theta)/((1+theta)(p+1)epsilon^(p+theta))]X^(p+2)
+O(X^(2p+2)).
```

For `p=4`, neither function produces a constant, linear physical term, or a
quadratic correction to the ADM Einstein-Hilbert operator. At `K=0`,
`Q_epsilon=1` and `V_epsilon=0`; the cosmological action reduces exactly to
the Einstein-Hilbert ADM action, up to the discarded boundary representative.

## 6. Local Matching Rule And PPN Prediction

The cosmological EFT is not extrapolated without a scale rule into resolved
Solar-System geometry. Define

```math
W=sqrt(abs(C_abcd C^abcd))/(H0/c)^2.
```

Let `B(x)=0` for `x<=0` and `B(x)=exp(-1/x)` for `x>0`. With
`W1=1e8`, `W2=1e9`, set

```math
s(W)=1                                  for W<=W1,
s(W)=B(1-t)/[B(1-t)+B(t)]               for W1<W<W2,
s(W)=0                                  for W>=W2,
t=(W-W1)/(W2-W1).
```

The complete low-energy action is

```math
S_eff=S_GR+integral d^4x s(W)[L_cos-L_GR].
```

It is a Wilsonian matching rule: in the `s=0` domain the action and all of its
variations are exactly those of GR with universal matter coupling. A flat FLRW
background has Weyl tensor zero and hence `s=1`. For a Schwarzschild Sun,

```math
W(r)=sqrt(48)(G M_sun/c^2)/[r^3(H0/c)^2].
```

The computed value at one AU is about `6e22`, far inside the GR domain. The
PPN prediction there is consequently

```math
gamma=beta=1, alpha1=alpha2=0.
```

The switch is deliberately C-infinity, so it has no artificial discontinuity at
the matching scale. As an EFT it is used only below its matching cutoff; no
claim about a UV completion is implied.

## 7. Tensor Sector

On the cosmological background the tensor quadratic action is

```math
S_T^(2)=(M_Pl^2/8) integral dt d^3x a^3 Q_epsilon(E)
[h_dot_ij^2-(partial_k h_ij)^2/a^2].
```

Therefore

```math
c_T^2=1,
Q_T(a)=Q_epsilon(E(a))>0,
dL_GW/dL_EM=sqrt(Q_T(0)/Q_T(z)).
```

The local matching does not alter the cosmological result because `W=0` on
FLRW.

## 8. Matter And CMB Interface

The complete input for linear matter perturbations is the action, not an
assumed `mu(a,k)` or a quasi-static approximation. The raw ADM derivatives
needed by a full 3+1 EFT implementation are obtained directly from

```math
L=Q_epsilon(X)[R3+K_ij K^ij-K^2]+2H0^2 V_epsilon(X).
```

In particular,

```math
L_R3=Q_epsilon,
L_Kij_R3=(Q_epsilon,X/(3H0))gamma^ij,
```

and

```math
L_Kij=2Q_epsilon(K^ij-K gamma^ij)
+[Q_epsilon,X(R3+K_mnK^mn-K^2)/(3H0)
+2H0 V_epsilon,X/3]gamma^ij.
```

Differentiating this last expression once more gives the full kinetic Hessian.
The generated `eft_coefficients.csv` records the background functions and
derivatives needed to implement these formulas without numerical
differentiation of a singular branch.

The Boltzmann backend must evolve the exact linear constraint system together
with photons, baryons, CDM, massive neutrinos, recombination, and lensing. The
protocol in `likelihood_pipeline.md` specifies the likelihood and rejection
rules.
