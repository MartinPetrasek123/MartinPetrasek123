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
3. Evaluate the exact ADM derivatives from `rfg_regularized.py` and write the
   table used by the 3+1 Einstein-Boltzmann module.
4. Evolve the full scalar, vector, and tensor linear system without a
   quasi-static approximation. The scalar constraints must be solved with the
   standard photon-baryon, neutrino, baryon, and CDM hierarchy present.
5. Reject a point when the module reports a ghost, gradient instability,
   singular constraint matrix, or non-positive tensor kinetic coefficient.
6. Export `C_ell^TT`, `C_ell^TE`, `C_ell^EE`, `C_L^phiphi`, `P(k,z)`,
   `f sigma8(z)`, and the standard-siren distance ratio.

H-EFTCAMB is an appropriate public backend because it evolves full linear
modified-gravity perturbations and exposes a full 3+1 EFT implementation with
stability checks. The RFG-R action derivatives, rather than a fitted CPL
surrogate, are the source of its cosmological functions.

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

The repository contains every RFG-R function and every deterministic rejection
condition. A full CMB likelihood requires two external public inputs that are
not vendored here: a compiled H-EFTCAMB installation and the official Planck
likelihood data package. The input/output contract is explicit:

```text
input:  generated/tables/eft_coefficients.csv
theory: full 3+1 RFG-R action derivative module
output: C_ell, P(k,z), f sigma8, and a stability flag
data:   Planck Legacy Archive likelihood packages
sampler: Cobaya nested or MCMC sampler
```

This boundary is operational, not theoretical: the action, variables,
parameters, output observables, datasets, and rejection rules are all fixed.
No empirical result is claimed until the external likelihood is actually run.

## 6. References For The Pipeline

1. B. Hu et al., Effective Field Theory of Cosmic Acceleration: an
   implementation in CAMB, Phys. Rev. D 89, 103530 (2014), arXiv:1312.5742.
2. G. Ye et al., H-EFTCAMB: A Cobaya-Integrated, Python-Wrapped Extension of
   EFTCAMB for Covariant Horndeski Gravity, arXiv:2603.01662 (2026).
3. Planck Collaboration V, Planck 2018 results. V. CMB power spectra and
   likelihoods, Astron. Astrophys. 641, A5 (2020), arXiv:1907.12875.
4. B. Bertotti, L. Iess, and P. Tortora, A test of general relativity using
   radio links with the Cassini spacecraft, Nature 425, 374-376 (2003).

