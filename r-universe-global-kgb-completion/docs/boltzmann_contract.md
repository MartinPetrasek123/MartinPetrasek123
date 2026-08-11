# Boltzmann implementation contract

The package supplies the complete background and Horndeski functions needed by
a linear Einstein--Boltzmann solver.  In the Bellini--Sawicki convention the
runtime functions are

\[
H(a)=H_0E(a),\qquad
\alpha_B(a),\qquad \alpha_K(a),\qquad
\alpha_M(a)=0,\qquad \alpha_T(a)=0,
\]

with standard minimally coupled baryons, CDM, photons, and neutrinos.  The
table `generated/ru_kgb_trajectory.csv` gives the exact values; the defining
functions are in `scripts/ru_kgb.py`, avoiding an interpolation ambiguity.

The native implementation is `scripts/generate_heftcamb_rph.py`, which maps
the Bellini--Sawicki braiding into the H-EFTCAMB RPH convention
`alpha_B^RPH=-alpha_B^BS/2`. `scripts/run_heftcamb_rph.py` executes the
coupled photon--baryon hierarchy, CDM, three massless standard neutrinos,
metric constraints, and KGB scalar with recombination, drag, lensed CMB
spectra, lensing-potential spectra, and linear matter power. It records the
201/401/601-node spline-convergence test in
`generated/heftcamb/convergence/heftcamb_convergence.json`.

The Planck wrapper `scripts/evaluate_planck_2018_fixed.py` evaluates raw
`C_ell` spectra against the official Planck 2018 Plik-lite TTTEEE, Commander
low-ell TT, SimAll low-ell EE, and CMB-dependent lensing likelihood files.
The reference calculation, including spline and perturbation-turn-on tests,
is collected in `generated/planck_2018_fixed_summary.json`. It is one fixed
parameter point with `A_planck=1`, not a sampled posterior or a joint
CMB--BAO--RSD likelihood.

With an H-EFTCAMB executable, an official Planck 2018 likelihood distribution,
and a Python environment containing `clipy` plus an HDF/FITS reader, the
external stages are reproduced as

```sh
HEFTCAMB_BIN=/path/to/camb python3 scripts/run_heftcamb_rph.py
CLIPY_SOURCE=/path/to/clipy PLANCK_2018_BASE=/path/to/plc_3.0 \
  /path/to/python-with-clipy scripts/evaluate_planck_2018_fixed.py \
  --spectra-dir generated/heftcamb/convergence/nodes_601 \
  --output generated/planck_2018_fixed_loglike_nodes_601.json
python3 scripts/summarize_planck_fixed_runs.py
```

The Planck data files are not redistributed by this package; the output report
records their local base directory and all spectrum conventions used in the
evaluation.

## Executed multi-redshift matter audit

`generate_heftcamb_rph.py --transfer-redshifts ...` writes the requested
transfer and linear total-matter power outputs in descending redshift order,
which is the forward-time order consumed by H-EFTCAMB. The exact selected KGB
point `alpha=0.0975`, `Omega_m0=0.3075` was executed at the central RSD
redshifts and at both sides of two finite-difference steps. The evaluator
`scripts/evaluate_kgb_rsd_native.py` computes

\[
\sigma_8^2(z)=\int d\ln k\,\frac{k^3P(k,z)}{2\pi^2}W^2(8k),
\qquad f\sigma_8=\frac{d\sigma_8}{d\ln a},
\]

directly from those H-EFTCAMB outputs and cross-checks the result against the
solver's velocity-density value `sigma8^2_vd/sigma8`. The produced records are
`generated/rsd_native_best_point/rsd_native_dz_0p005.json` and
`generated/rsd_native_best_point_dz_0p0025/rsd_native_dz_0p0025.json`.
They are a numerical and physical audit, not an RSD likelihood: a full survey
calculation requires the survey covariance, window/AP mapping, and nonlinear
nuisance model and is not represented by the compact local data file.
