# Exact RFG-R ADM-to-Extended-EFT Map

This note maps the RFG-R cosmological action, not a fitted background history,
to the extended EFT basis used for linear perturbations. The reference mapping
is Frusciante, Papadomanolakis and Silvestri, arXiv:1601.04064, Eqs. (7)-(19).
The reference uses the opposite extrinsic-curvature sign, so this map makes
the explicit conversion `K_reference=-K_RFG`; RFG-R itself has
`K_ij=H gamma_ij` on FLRW.

## ADM input

With `S=K_ij K^ij`, `R3={}^{(3)}R`, and `X=-K_reference/(3 H0)`, the gravitational ADM
Lagrangian is

```math
L = M_Pl^2/2 [Q(X)(R3 + S - K^2) + 2 H0^2 V(X)].
```

It has no explicit lapse dependence, so `L_N=L_NN=L_NR=L_NS=0`. On a flat
FLRW background, `R3=0`, `S=3H^2`, and `K=3H`. The nonzero derivatives needed
by the general ADM map are

```math
E = L_R = M_Pl^2 Q/2,
C = L_KR = -M_Pl^2 Q_X/(6 H0),
L_S = M_Pl^2 Q/2,
```

and

```math
A = M_Pl^2/2[-2 X^2 Q_XX/3 - 8 X Q_X/3 - 2Q + 2V_XX/9].
```

The remaining background terms are evaluated directly from the reconstructed
potential and the implicit RFG-R branch.

## Result

The exact mapped functions include

```math
Omega = Q-1,
bar_M1^3 = -M_Pl^2 dot(Q),
bar_M2^2/M_Pl^2 = X^2 Q_XX/3 + 4 X Q_X/3 - V_XX/9,
bar_M3^2 = 0,
bar_m5 = -M_Pl^2 Q_X/(3H0).
```

The last coefficient multiplies the extended action operator

```math
bar_m5/2 delta R3 delta K = -M_Pl^2 Q_X/(6H0) delta R3 delta K.
```

It is forced by the `Q(K) R3` term in the RFG-R action. It cannot be set to
zero without changing the theory. The table records the basis-level
coefficients above rather than solver-specific gamma normalizations. The usual
six-function EFTCAMB input, whatever convention is used for its gamma
functions, does not encode `bar_m5 delta R3 delta K`.

## Executable calculation

Run

```bash
python3 scripts/validate_extended_eft_mapping.py
python3 scripts/extended_eft_mapping.py
```

The machine-readable result is `generated/tables/extended_eft_mapping.csv`.
Values use the dimensionless normalizations stated in the script header. The
validation independently finite-differences both implicit-background
derivatives and checks the exact `bar_m5=-M_Pl^2 Q_X/(3H0)` identity.

## Backend consequence

The downloaded and compiled public H-EFTCAMB revision exports only the usual
`Omega`, `gamma_1` through `gamma_6` operator set. A source search confirms
that it has no `bar_m5 delta R3 delta K` coefficient. Therefore it is suitable
for an unmodified GR reference calculation but cannot yield an exact RFG-R
CMB, lensing, or matter prediction. An RFG-R Boltzmann implementation must
add the extended operator to the scalar equations and its stability checks,
then validate a GR limit before an observational likelihood is run.
