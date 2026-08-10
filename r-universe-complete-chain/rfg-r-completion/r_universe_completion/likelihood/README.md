# External Likelihood Runtime

This directory intentionally does not vendor either Planck data or a modified
Einstein-Boltzmann code. Both are third-party distributions with their own
licenses and release procedures.

## Required Runtime

1. Clone H-EFTCAMB with submodules and build its Python interface:

```bash
git clone --recursive https://github.com/EFTCAMB/EFTCAMB.git
cd EFTCAMB/fortran
make python
```

The build requires a Fortran compiler, C compiler, BLAS/LAPACK, NumPy, SciPy,
SymPy, and packaging. The public H-EFTCAMB README documents the platform-specific
requirements.

2. Reproduce the pure-gravity gate in
`../scripts/validate_extended_eft_scalar_stability.py`. It finds a degenerate
standalone scalar kinetic coefficient, so a one-scalar backend is not an
RFG-R implementation.

3. Implement the already derived sourced photon--baryon--CDM--neutrino RFG-R
   system, including `bar_m5 delta R3 delta K`, from
   `../docs/photon_baryon_cdm_neutrino_reduction.md` and
   `../scripts/multifluid_reduction.py`. Its finite constraint matrix, exact
   kinetic hierarchy, initial conditions, and GR limit must be preserved in
   the solver before spectra are computed.

4. Install the official Planck 2018 likelihood code and data from the Planck
Legacy Archive, then connect them to Cobaya.

## Mandatory Per-Sample Gate

The RFG-R module must return a failure before a likelihood is evaluated when
the full multi-fluid scalar system has a ghost, a gradient instability, a
singular constraint matrix, or `Q_T<=0` on the integration domain. A backend
that replaces the degenerate pure-gravity result with a guessed kinetic
operator, or only uses the background table, a CPL fit, a quasi-static
approximation, or a compressed CMB distance prior is not the RFG-R likelihood
defined by this package.

## Output Contract

```text
input:  extended_eft_mapping.csv, eft_coefficients.csv, and sampled cosmological parameters
output: C_ell TT/TE/EE, C_L phiphi, P(k,z), f sigma8(z), dL_GW/dL_EM(z)
gate:   full multi-fluid scalar and tensor stability report
data:   official Planck likelihood plus selected BAO/SN/RSD/GW data
sampler: Cobaya
```
