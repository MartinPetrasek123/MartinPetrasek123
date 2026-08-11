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

## Conditional Planck profile

The package now includes a declared coarse-to-local grid in the two action
parameters `alpha` and `Omega_m0`, using native H-EFTCAMB spectra and the
official Planck 2018 Plik-lite TTTEEE, Commander, SimAll, and lensing
likelihood objects.  The final local 601-node grid has its minimum at
`alpha=0.04`, `Omega_m0=0.3085`, with
`-2 ln L = 1075.3470998357473` at fixed `A_planck=1`.  Its independent action
validation is recorded in
`generated/planck_profile_final_local/best_point_action_validation.json`.

For the same `Omega_m0`, radiation, primordial inputs, and fixed spectra, the
standard LCDM reference has `-2 ln L = 1066.1743939041603` at `A_planck=1`.
Profiling only the Planck absolute-calibration nuisance on the declared
`A_planck` grid gives minima `1066.1369110626738` for KGB and
`1055.786875917371` for LCDM, both at `A_planck=1.0025`.  Thus the local
conditional KGB point is higher by `10.350035145302854` in `-2 ln L` than the
matched fixed LCDM point.  This is a diagnostic comparison at fixed remaining
parameters, not a posterior, evidence calculation, or model-selection claim.

## Boundary of the conclusion

These computations establish a globally specified, classically regular
covariant scalar-tensor candidate on the stated branch. A native
Einstein--Boltzmann evolution, an official Planck 2018 fixed-point CMB plus
lensing likelihood, and a limited conditional two-parameter grid have been
executed. They do **not** derive this KGB action from the primitive RCD
capacity field, provide a Planck posterior, optimize the full cosmological and
nuisance parameter space, perform a non-linear screening solution, or supply
a quantum UV completion. Those remaining statements are not implied by this
document.
