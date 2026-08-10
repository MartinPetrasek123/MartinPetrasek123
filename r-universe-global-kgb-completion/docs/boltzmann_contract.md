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
