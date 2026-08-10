# Physical status of the global KGB action

## What this package defines

The action

\[
S=\int d^4x\sqrt{-g}\left[
\frac{M_{\rm Pl}^2}{2}R+A(\phi)X+B(\phi)X^2-V(\phi)-C(\phi)X\Box\phi
\right]+S_m[g,\psi_m]
\]

is a covariant `G2+G3` Horndeski/kinetic-gravity-braiding action. For this
action, diffeomorphism invariance gives the Noether identity relating the
metric Euler equation and the scalar Euler equation.

As a `G2+G3` Horndeski subclass, the action has the standard two tensor modes
and one scalar mode, with no higher-derivative Ostrogradsky mode.  Positivity
of the tensor kinetic term, scalar kinetic coefficient, and scalar gradient
coefficient is still a branch condition and is checked numerically here.

## Executed existence gates

On `1e-8 <= a <= 1e3`, the reproducible validation obtains positive
`N_s`, `D`, `Q_s`, `F_s`, and luminal positive tensor coefficients.  The
independent action-identity validator additionally checks the field-coordinate
background derivative, R-sector derivative, `alpha_B`, `alpha_K`, and the
homogeneous scalar equation using finite differences and the covariant action
formulae, rather than the reconstruction residuals alone.

The future saturation is part of the theory definition: `a_saturation=2` and
the target scalar sound speed `c_s^2=1` are prescribed action choices. They
are not observational inferences. The regular early-time braiding closure
`alpha_B=Omega_R/3` is fixed by the radiation-limit scalar normalization; it
is not a sampled cosmological parameter or a derivation from the RCD capacity
field.

## Boundary of the conclusion

These computations establish a globally specified, classically regular
covariant scalar-tensor candidate on the stated branch. A native
Einstein--Boltzmann evolution and an official Planck 2018 fixed-point CMB plus
lensing likelihood have been executed. They do **not** derive this KGB action
from the primitive RCD capacity field, provide a Planck posterior or model
comparison, perform a non-linear screening solution, or supply a quantum UV
completion. Those remaining statements are not implied by this document.
