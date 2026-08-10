# R-Universe RFG-R Completion

This package defines RFG-R: a regular, multiscale effective completion of the
R-Universe preferred-foliation branch. It is designed to answer a precise
question: what theory is to be tested when one asks for matter growth, CMB
spectra, and Solar-System consistency?

RFG-R is not a rhetorical strengthening of the original branch. It specifies:

- one smooth nonlinear action on the cosmological branch;
- a reconstructed potential that preserves the R-Universe background to
  `O((epsilon/X)^p)`;
- a local GR matching domain that makes the PPN prediction unambiguous;
- universal minimal coupling of matter to the Jordan metric;
- an exact sourced lapse--shift reduction for one canonical scalar field,
  including a General-Relativity normalization check;
- a pure-gravity extended-EFT scalar audit, which finds a degenerate scalar
  kinetic coefficient and therefore blocks a direct one-scalar Boltzmann run;
- an exact photon--baryon--CDM--neutrino finite constraint reduction and an
  untruncated photon-polarization/neutrino kinetic hierarchy;
- a full linear Einstein-Boltzmann likelihood protocol, still awaiting a
  solver implementation and an executed data run;
- an exact ADM-to-extended-EFT map, including the nonzero
  `bar_m5 delta R3 delta K` coefficient required by the action.

The resulting physical picture is simple: a relational response is active on
homogeneous Hubble-scale geometry, while the locally resolved weak-field EFT is
General Relativity. The crossover is explicit rather than implicit.

## Model Definition

For `p=4`, `0<theta<2`, and `epsilon>0`, define

```math
R_epsilon(X)=Omega_R0 X^(p+2)/(X^p+epsilon^p)^(1+theta/p),
```

```math
Q_epsilon(X)=1-[Omega_R0/(1+theta)]X^p/(X^p+epsilon^p)^(1+theta/p).
```

The potential is fixed, not guessed:

```math
V_epsilon(X)=-3X integral_0^X ds/s^2 [s^2-R_epsilon(s)-s^2Q_epsilon(s)-s^3Q_epsilon'(s)].
```

The cosmological action is the RFG action with `Q,V` replaced by these
functions. The background equation is

```math
E^2=Omega_m0 a^-3+Omega_r0 a^-4+R_epsilon(E).
```

The local matching action is exactly Einstein-Hilbert plus minimally coupled
matter whenever the dimensionless Weyl indicator

```math
W=sqrt(abs(C_abcd C^abcd))/(H0/c)^2
```

exceeds `W_LOCAL_GR=1e9`. A C-infinity switching function is given in
`docs/completion_derivation.md`. FLRW has `W=0`; at one astronomical unit from
the Sun the code obtains `W about 6e22`.

## Contents

- `paper/R_Universe_RFG_R_Completion.pdf` - standalone paper.
- `docs/completion_derivation.md` - full mathematical definition and limits.
- `docs/likelihood_pipeline.md` - matter, CMB, and PPN likelihood protocol.
- `docs/extended_eft_mapping.md` - exact action-to-extended-EFT map and the
  public-backend compatibility boundary.
- `docs/extended_eft_scalar_audit.md` - full pure-gravity scalar-action audit
  retaining `bar_m5`, including the numerical degeneracy result.
- `docs/photon_baryon_cdm_neutrino_reduction.md` - exact finite
  photon--baryon--CDM--neutrino constraint reduction plus the untruncated
  photon-polarization and neutrino kinetic hierarchy.
- `docs/canonical_scalar_reduction.md` - exact canonical-scalar quadratic
  reduction and its executable stability diagnostic.
- `scripts/rfg_regularized.py` - background solver, potential reconstruction,
  and ADM coefficient table.
- `scripts/validate_completion.py` - independent numerical checks.
- `scripts/extended_eft_mapping.py` - executable extended-EFT coefficient
  table, including `bar_m5`.
- `scripts/validate_extended_eft_mapping.py` - independent mapping checks.
- `scripts/extended_eft_scalar_stability.py` - complete pure-gravity extended
  EFT scalar audit; deliberately returns failure on a degenerate action.
- `scripts/validate_extended_eft_scalar_stability.py` - regression test for
  that physical failure result.
- `scripts/multifluid_reduction.py` - exact Sorkin--Schutz constraint blocks,
  Schur reduction, Thomson terms, massless hierarchy, and massive-neutrino
  phase-space interface.
- `scripts/validate_multifluid_reduction.py` - rational GR rank test and the
  425-point RFG-R multi-fluid core audit.
- `scripts/ppn_likelihood.py` - local GR matching and Cassini likelihood.
- `generated/` - tables and figures generated from the scripts.

## Reproduce

```bash
bash scripts/run_all.sh
```

The command checks the exact local limit, high-X recovery, potential
reconstruction, background closure, positive tensor normalization, and the PPN
domain. It also checks the canonical-scalar General-Relativity limit and its
reference sourced-constraint diagnostic, reproduces the pure-gravity
extended-EFT scalar degeneracy, and runs the exact multi-fluid core audit.

## Empirical Rule

The included CMB/matter protocol is a full-spectrum likelihood definition, not
a compressed-distance surrogate. The exact action map has a nonzero extended
`bar_m5 delta R3 delta K` coefficient. Its finite physical multi-fluid
constraint reduction and its untruncated kinetic hierarchy are now supplied,
but the public stock H-EFTCAMB interface still does not expose `bar_m5`; it is
only a GR reference calculation, not an RFG-R prediction engine. The package
never substitutes a background fit for a CMB fit or claims empirical
preference over LCDM.
