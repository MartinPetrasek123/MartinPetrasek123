# Matter, CMB, And PPN Likelihood Pipeline

## 1. Claim Tested

The pipeline tests RFG-R as a specified low-energy theory against the same
observables used for LCDM. It does not infer a CMB conclusion from `H(z)` or
from a distance-prior compression.

The sampled cosmological parameters are

```text
omega_b, omega_cdm, H0, ln(10^10 A_s), n_s, tau_reio,
sum_mnu, theta, epsilon, p.
```

For the published reference branch `p=4` and `epsilon=1e-8` are fixed because
they affect neither CMB-scale nor background observables at machine precision.
They remain explicit model parameters so the local completion is reproducible.

## 2. Theory Evaluation

For each proposed point:

1. Solve `E^2=Omega_m0 a^-3+Omega_r0 a^-4+R_epsilon(E)` on the positive branch.
2. Reject if a root is absent, density closure fails, or `Q_T<=0` anywhere on
   `a in [1e-8,1]`.
3. Evaluate the exact ADM derivatives from `rfg_regularized.py`, map them
   with `extended_eft_mapping.py`, and write
   `generated/tables/extended_eft_mapping.csv`. The map contains the nonzero
   extended coefficient `bar_m5=-M_Pl^2 Q_X/(3H0)`.
4. Reduce the complete extended EFT *pure-gravity* scalar action, including
   ``bar_m5 delta R3 delta K``, then use the exact sourced finite reduction in
   `docs/photon_baryon_cdm_neutrino_reduction.md`. Its lapse--shift equations
   must be solved with the density and momentum moments of every hierarchy.
5. Evolve the full scalar, vector, and tensor system without a quasi-static
   approximation. For the scalar sector this means the exact photon
   temperature/polarization hierarchy and collisionless neutrino hierarchy;
   no fluid closure or compressed CMB observable is permitted.
6. Reject a point when the module reports a ghost, gradient instability,
   singular constraint matrix, a zero of the derived curvature coefficient
   `mu_zeta=M_00-dot(B_00)`, or non-positive tensor kinetic coefficient.
7. Export `C_ell^TT`, `C_ell^TE`, `C_ell^EE`, `C_L^phiphi`, `P(k,z)`,
   `f sigma8(z)`, and the standard-siren distance ratio.

The audited public H-EFTCAMB source is useful for checking the standard
six-function EFT interface. It is not an exact RFG-R backend: its exposed
`Omega`, `gamma_1` through `gamma_6` input set has
no `bar_m5 delta R3 delta K` operator. Setting that coefficient to zero would
change the RFG-R action. The RFG-R Boltzmann module must therefore add this
extended operator to the scalar equations and stability checks before any
RFG-R spectrum or likelihood is evaluated. The action derivatives, rather
than a fitted CPL surrogate, remain the source of its cosmological functions.

## 2.1 Executed RFG-R Gate

The exact pure-gravity reduced action was evaluated by
`scripts/extended_eft_scalar_stability.py` on
`a in [10^-7,1]` and `k/H0 in [10^-4,10^5]`.  It retains the mapped nonzero
`bar_m5` coefficient and the nonzero `W7` term.  The lapse/shift constraint
discriminant is positive, but the coefficient of `zeta_dot^2` cancels on all
2,401 sampled points.  Hence the *pure-gravity* scalar quadratic action is
degenerate and its standalone scalar sound speed is undefined.  This is a
strong-coupling/closure gate, not a successful stability result.

The exact multi-fluid reduction now identifies the null direction as the
algebraic curvature constraint after
`Delta_i=delta_i+3(1+w_i)zeta`.  Its coefficient
`mu_zeta=M_00-dot(B_00)` is derived with analytic action derivatives in
`scripts/rfg_dae_closure.py`.  It crosses zero on the reference branch: the
generated 49-point table has 24 root-bearing scale factors and gives
`k/H0=2.51545672221` at `a=1`.  The zero is a singularity of the exact linear
DAE, not an unimplemented collision or hierarchy term.  Therefore the
reference RFG-R action must be rejected before spectra, a sampler, or any
likelihood component is evaluated.  Adding a gravitational kinetic operator
to avoid it would define a different theory and would have to be stated,
derived, and tested as such.

## 2.2 RFG-RXi Decision Gate

RFG-RXi is not a numerical repair of the rejected RFG-R branch.  Its action
adds the separately derived background-null operator

```text
M_Pl^2/2 int d^4x N sqrt(gamma) Xi [R3 + sigma_ij sigma^ij].
```

Xi is a new dimensionless EFT coefficient.  The package does not infer it
from data: the unfitted checks at Xi=1 and Xi=2 are two declared benchmark
evaluations.  Both retain the R-Universe background exactly and return no
sampled `mu_zeta` root on their separate 49-by-81 `(a,k)` audits.  This is a
necessary internal gate, not a spectrum, posterior, or likelihood value.

An RFG-RXi Boltzmann module must start from its own augmented ADM derivatives,
repeat the lapse--shift and spatial-metric reduction, and use the resulting
Xi-dependent hierarchy throughout recombination and propagation.  It may not
reuse an RFG-R transfer function across the rejected curvature-constraint
surface or declare the GR/CAMB reference to be an RFG-RXi prediction.

### Directly Evaluated RFG-RXi Factors

The local Weyl matching switch is exactly zero in the Solar-System domain.
The RFG-RXi operator is therefore absent from the local action, giving
`gamma=beta=1`, `alpha1=alpha2=0`, and the directly evaluated Cassini factor
`-2 ln L_Cassini=0.8336483932` for both Xi=1 and Xi=2. The same action gives
the exact tensor relation

```text
Q_T^Xi(a)=Q_epsilon(a)+Xi,
c_T=1,
dL_GW/dL_EM=sqrt(Q_T^Xi(1)/Q_T^Xi(a)).
```

At the stated Planck-recorded background input, the deterministic prediction
at z approximately 1 is `0.9565968846` for Xi=1 and `0.9717969466` for Xi=2.
These are action-defined partial factors, not a combined cosmological
likelihood or a comparison with LambdaCDM. The generator, validation, and
tables are `scripts/rfg_xi_observables.py`,
`scripts/validate_rfg_xi_observables.py`,
`generated/tables/rfg_xi_observables.csv`, and
`generated/tables/rfg_xi_validated_factors.csv`.

## 2.3 Executed GR Infrastructure Reference

The standard photon/baryon/CDM/neutrino, recombination, lensing and transfer
machinery was executed with the pinned public `camb==2.0.1` wheel at the
published Planck 2018 TT,TE,EE+lowE+lensing posterior means. The generated
reference has `sigma8(z=0)=0.8110325278646`; lensed CMB, lensing-potential and
linear-matter samples are written to `generated/tables/gr_camb_reference.csv`.
`scripts/validate_gr_reference_camb.py` reproduces the pinned values.

Official Planck low-T, low-E and lensing components were also executed through
Cobaya 3.6.2 using the Planck BBN-consistency preset and the defined nominal
calibration `A_planck=1`. The individual values and the combined
`-2 ln L=428.3415086187` are recorded in
`generated/tables/gr_planck_2018_lowell_lensing.csv`; their regression script
is `scripts/validate_gr_planck_lowell_lensing.py`.

This is a GR data-interface test only. Stock CAMB has `bar_m5=0`, whereas
RFG-R has a nonzero mapped coefficient. It does not produce an RFG-R spectrum,
posterior, best fit or comparison to LCDM. The official high-ell Plik package
is installed but not assigned guessed calibration/foreground values; it must
enter an actual joint fit. Details and commands are in
`docs/gr_camb_reference.md`, `docs/planck_data_reference.md`, and
`docs/backend_capability_audit.md`.

## 3. Likelihoods

The baseline joint log likelihood is

```math
ln L_total=ln L_Planck18+ln L_Planck_lensing+ln L_BAO
+ln L_SN+ln L_RSD+ln L_PPN+ln L_GW.
```

Required baseline components:

- Planck 2018 low-ell temperature and polarization, high-ell TT/TE/EE, and
  lensing likelihood from the official Planck Legacy Archive.
- A BAO likelihood selected without using LCDM-compressed observables outside
  their stated validity range.
- Supernova and RSD likelihoods that consume the direct theory predictions.
- The PPN likelihood below.
- Standard-siren events with electromagnetic redshift, using
  `dL_GW/dL_EM=sqrt(Q_T(0)/Q_T(z))` in the cosmological domain.

## 4. PPN Likelihood

In the matched local GR domain the prediction is

```math
gamma=beta=1, alpha1=alpha2=0.
```

The Cassini factor implemented in `scripts/ppn_likelihood.py` is

```math
-2 ln L_Cassini=[(gamma-1-2.1e-5)/(2.3e-5)]^2.
```

The reference RFG-R point gives `-2 ln L_Cassini=0.833648...`. Additional PPN
measurements can be added as Gaussian factors with their published covariance.

## 5. Reproducible Execution Boundary

The repository contains every RFG-R background function, the exact
ADM-to-extended-EFT map, the finite physical constraint reduction, the
untruncated kinetic hierarchy, and the action-derived DAE closure.  The last
step rejects the published reference branch because its curvature constraint
loses rank at finite wavenumber.  The following input/output contract applies
only to a newly specified action that first passes the same `mu_zeta` gate:

```text
input:  generated/tables/extended_eft_mapping.csv
theory: full 3+1 RFG-R action-derivative module extended by bar_m5
output: C_ell, P(k,z), f sigma8, and a stability flag
data:   Planck Legacy Archive likelihood packages
sampler: Cobaya nested or MCMC sampler
```

The boundary is therefore physical as well as numerical: no RFG-R spectra,
posterior, or CMB/matter likelihood is valid for the published reference
branch. The same data-interface boundary applies to RFG-RXi until its own
action-faithful solver and likelihood have been run. The executed GR reference
does not change either statement.

## 6. References For The Pipeline

1. B. Hu et al., Effective Field Theory of Cosmic Acceleration: an
   implementation in CAMB, Phys. Rev. D 89, 103530 (2014), arXiv:1312.5742.
2. N. Frusciante, G. Papadomanolakis and A. Silvestri, An Extended action for
   the effective field theory of dark energy: stability analysis and a complete
   guide to the mapping at the basis of EFTCAMB, arXiv:1601.04064 (2016).
3. G. Ye et al., H-EFTCAMB: A Cobaya-Integrated, Python-Wrapped Extension of
   EFTCAMB for Covariant Horndeski Gravity, arXiv:2603.01662 (2026). This is
   a reference-backend citation, not evidence that its public operator basis
   contains the RFG-R extended coefficient.
4. Planck Collaboration V, Planck 2018 results. V. CMB power spectra and
   likelihoods, Astron. Astrophys. 641, A5 (2020), arXiv:1907.12875.
5. B. Bertotti, L. Iess, and P. Tortora, A test of general relativity using
   radio links with the Cassini spacecraft, Nature 425, 374-376 (2003).
6. A. De Felice, N. Frusciante, and G. Papadomanolakis, On the stability of
   Horndeski theories with matter fields, JCAP 03, 027 (2017),
   arXiv:1609.03599.
7. J. Hwang and H. Noh, Cosmic microwave background anisotropies and
   power spectra from a unified gauge-ready formulation, Phys. Rev. D 65,
   023512 (2001), arXiv:astro-ph/0102005.
