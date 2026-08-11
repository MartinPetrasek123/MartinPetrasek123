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

## Executed conditional CMB and late-time record

The package contains a 154-point fixed-input surface built from native
H-EFTCAMB spectra and the official Planck 2018 Plik-lite TTTEEE, Commander,
SimAll, and lensing likelihood objects, added to exact-background Pantheon+
full-covariance, DESI DR2 BAO full-covariance, and chronometer evaluations.
At fixed `A_planck=1` its lowest executed conditional ordinate is
`2493.222580514448` at `alpha=0.1`, `Omega_m0=0.3075`; this is a grid result,
not a continuous optimum.

A separate refinement selects `alpha=0.0975`, `Omega_m0=0.3075` and validates
the reconstructed action at that point. Its Planck calibration grids have
minima `1063.6712514661913` for KGB at `A_planck=1.00275` and
`1060.2194088862268` for matched LCDM at `A_planck=1.00325`. With the
solver-derived `r_drag=147.32 Mpc`, the independently integrated late-time
statistics are `1416.5609849186849` (KGB) and `1422.5765927552` (LCDM).
The resulting executed conditional sums are `2480.232236384876` and
`2482.7960016414268`, respectively: KGB minus LCDM is `-2.563765256550596`.
Every non-displayed cosmological, primordial, recombination, and high-ell
nuisance input is fixed, so this is neither a posterior nor evidence. Its
machine-readable ledger is `generated/kgb_multprobe_conditional_summary.json`.

The selected KGB point also passes the same analytic cubic-KGB local screening
gate with `r_V=102.55866003368482 pc` and a conservative
`|gamma_PPN-1|=1.0462269740540104e-15` envelope. This is not an ephemeris fit.

## Native RSD audit

The selected KGB point was evolved natively at 58 transfer redshifts. From its
linear total-matter `P(k,z)`, `sigma8` is directly integrated and `f sigma8`
is obtained by symmetric `d sigma8/d ln a`. It covers all 23 rows of the local
RSD compilation. The independent H-EFTCAMB velocity-density output agrees with
the power-spectrum derivative to at most `9.385121718008538e-05`; halving the
redshift step changes the largest predicted `f sigma8` by
`5.00883760279347e-05`.

This is deliberately a residual audit, not an RSD likelihood: the supplied
WiggleZ and SDSS-IV groups lack their survey covariance, window functions,
fiducial-cosmology/AP mapping, and nonlinear nuisance treatment. Its
diagnostic-only diagonal residual sum (`18.947293535743064` for the 16 rows
labelled `diag`) is not combined with the conditional CMB--late-time sum.

## Boundary of the conclusion

These computations establish a globally specified, classically regular
covariant scalar-tensor candidate on the stated branch, with a native
Einstein--Boltzmann execution, an official Planck likelihood, an exact
late-time geometric conditional calculation, a local screening gate, and a
native RSD audit. They do **not** derive this KGB action from the primitive RCD
capacity field, provide a joint posterior/evidence, include a survey-complete
RSD or weak-lensing likelihood, optimize the full cosmological and nuisance
parameter space, perform a non-linear screening solution, or supply a quantum
UV completion. Those remaining statements are not implied by this document.
