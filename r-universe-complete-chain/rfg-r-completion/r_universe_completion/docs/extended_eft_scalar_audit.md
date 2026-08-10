# RFG-R Extended-EFT Pure-Gravity Scalar Audit

This note evaluates the scalar action implied by the exact RFG-R
ADM-to-extended-EFT map. It is not a matter, CMB, or likelihood calculation.
Its purpose is narrower and prior: determine whether the *pure-gravity*
extended-EFT scalar sector can be evolved as an ordinary standalone scalar.

## Input action and coefficients

The input is the exact map in `extended_eft_mapping.py`, including

```math
bar_m5 = -M_Pl^2 Q_X/(3 H0),
bar_M3^2 = 0,
m_2^2 = lambda_i = hat_M^2 = 0.
```

The nonzero `bar_m5 delta R3 delta K` term must be retained. For the
unitary-gauge scalar action of Frusciante, Papadomanolakis and Silvestri,
arXiv:1601.04064, Eqs. (85)-(86), the relevant coefficients are

```math
W0 = -[Q + 3 H bar_m5 + 3 dot(bar_m5)]/a^2,
W1 = c + 2 M2_4 - 3 H^2 Q - 3 H dot(Q)
     - 3 H^2 bar_M3^2/2 - 9 H^2 bar_M2^2/2 - 3 H bar_M1^3,
W4 = [-2 H Q-dot(Q)-H bar_M3^2-bar_M1^3-3H bar_M2^2]/a^2,
W5 = [2Q+bar_M3^2+3bar_M2^2]/a^2,
W6 = [-2Q-6H bar_m5]/a^2,
W7 = -[bar_M3^2+bar_M2^2]/(2a^4).
```

All values are evaluated with `H0=M_Pl=1` and the exact implicit RFG-R
background. The lapse--shift constraint discriminant is

```math
D = W4^2 - 4 W1 W7.
```

## Reduced kinetic coefficient

With the RFG-R values `m_2^2=lambda_i=0`, the coefficient of
`zeta_dot^2` after eliminating lapse and shift is

```math
L_zeta_dot_zeta =
 [(6a^2 W7+W5)(3a^4 W4^2+2a^2 W1 W5)]/[2a^2 D].
```

This is the `L_{dot zeta dot zeta}` specialization of Eq. (111) of the cited
extended-EFT derivation. A nonzero `bar_m5` remains present in the reduced
gradient and mixing terms; it cannot be omitted from the calculation merely
because it does not enter this kinetic expression.

## Executed result

Run:

```bash
python3 scripts/extended_eft_scalar_stability.py
python3 scripts/validate_extended_eft_scalar_stability.py
```

The first command emits the full grid to
[`extended_eft_scalar_stability.csv`](../generated/tables/extended_eft_scalar_stability.csv).
It deliberately exits nonzero because a zero kinetic coefficient is a physics
gate. The second command is the regression test that verifies the known result
and exits zero.

For the reference branch the audit uses 49 logarithmic scale factors over
`10^-7 <= a <= 1` and 49 logarithmic wavenumbers over
`10^-4 <= k/H0 <= 10^5`, for 2,401 points. It finds

```text
min constraint discriminant       = 2.514017848646e+00
kinetic-degenerate points        = 2401 / 2401
max relative kinetic numerator   = 4.357e-16
```

Thus the cancellation in the numerator of `L_zeta_dot_zeta` is numerical-zero
to double-precision accuracy on the full audited branch. The standalone
pure-gravity scalar sound speed is undefined; treating the ratio `G/L` as a
finite large number would be a numerical artefact.

## Consequence and limit of the audit

This result blocks a direct transfer of the RFG-R map to a conventional
one-scalar EFT/Boltzmann solver, in addition to the fact that the public
H-EFTCAMB interface has no `bar_m5 delta R3 delta K` input. It does **not** by
itself prove that the physical matter-coupled theory has no well-posed scalar
sector: matter sources the lapse and shift constraints. The package's separate
one-canonical-scalar calculation already illustrates that the sourced
reduction need not have the same kinetic rank.

The exact photon--baryon--CDM--neutrino finite constraint action and kinetic
hierarchy are now supplied in
[`photon_baryon_cdm_neutrino_reduction.md`](photon_baryon_cdm_neutrino_reduction.md),
including the mapped `bar_m5` operator. The remaining task is to implement
that differential--algebraic hierarchy, its recombination and initial
conditions, in a solver; then its integrated gradients, spectra and likelihood
can be evaluated. Adding a new gravitational kinetic operator instead would
define a new theory, not complete this RFG-R action.
