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
   singular constraint matrix, or non-positive tensor kinetic coefficient.
7. Export `C_ell^TT`, `C_ell^TE`, `C_ell^EE`, `C_L^phiphi`, `P(k,z)`,
   `f sigma8(z)`, and the standard-siren distance ratio.

The downloaded public H-EFTCAMB revision is useful as a compiled GR reference
and for checking the standard six-function EFT interface. It is not an exact
RFG-R backend: its exposed `Omega`, `gamma_1` through `gamma_6` input set has
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

Consequently the reference RFG-R branch cannot be passed directly to a
single-scalar EFT/Boltzmann backend. The physical
photon--baryon--CDM--neutrino reduction from the same action is now recorded in
`docs/photon_baryon_cdm_neutrino_reduction.md`: its finite monopole/dipole core
has four positive material directions and one constraint, while the photon and
neutrino anisotropic-stress sectors remain exact kinetic hierarchies. The
existing one-canonical-scalar diagnostic is not a substitute for that system.
Adding a gravitational kinetic operator would define a different theory and
would have to be stated, derived, and tested as such.

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
ADM-to-extended-EFT map, the finite physical constraint reduction, and the
untruncated kinetic hierarchy. The present action does not pass directly to a
one-scalar backend. A full CMB likelihood still requires an implementation of
this action-defined differential-algebraic hierarchy in a solver retaining
`bar_m5`, followed by the official Planck likelihood data package. The
input/output contract for that calculation is:

```text
input:  generated/tables/extended_eft_mapping.csv
theory: full 3+1 RFG-R action-derivative module extended by bar_m5
output: C_ell, P(k,z), f sigma8, and a stability flag
data:   Planck Legacy Archive likelihood packages
sampler: Cobaya nested or MCMC sampler
```

The remaining boundary is numerical/data-facing: no RFG-R spectra, posterior,
or CMB/matter likelihood has been executed. No empirical result is claimed.

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
