# R-Universe global R-alpha KGB completion

This package turns the R-alpha background already fitted in the accompanying
late-time R-Universe manuscript into one explicit, globally defined,
luminal KGB action.  It replaces the missing `kgb_reconstruction.py` named by
that manuscript with reproducible code and closes the scalar no-ghost and
no-gradient conditions analytically and numerically.

Run:

```sh
bash scripts/run_all.sh
```

Outputs:

- `generated/ru_kgb_trajectory.csv` - action and EFT functions;
- `generated/independent_action_identities.json` - independently recomputed action identities;
- `generated/validation.json` - all numerical gates;
- `generated/stability_scan.csv` and `.json` - physical stability scan around the fit;
- `generated/matter_qs_prediction.csv` - derived high-k linear matter response and growth;
- `generated/ppn_screening.json` - solar-system cubic-KGB screening gate;
- `generated/cmb_prerecombination_gate.json` - sound-horizon and pre-recombination EFT inputs;
- `generated/planck_2018_fixed_summary.json` - official Planck 2018 fixed-point CMB plus lensing likelihood and numerical checks;
- `generated/planck_profile_final_local/summary.json` - final local conditional KGB Planck grid with matched fixed LCDM references;
- `generated/planck_profile_final_local/best_point_action_validation.json` - action-level validation at the local grid minimum;
- `generated/planck_profile_final_local/calibration_kgb_refined/summary.json` and `calibration_lcdm_refined/summary.json` - fixed-spectrum Planck-calibration profiles;
- `generated/ru_kgb_stability.pdf` and `.png` - stability figure.

The standalone technical paper source is
`paper/R_Universe_Global_KGB_Completion.tex`. It references this exact source
tree and can be included as the covariant-completion companion to the
late-time R-Universe manuscript.

The construction preserves the R-alpha background exactly for all observed
redshifts. It includes a derived subhorizon matter prediction, a
decoupling-limit PPN screening gate, and a native H-EFTCAMB calculation of the
photon--baryon--CDM--massless-neutrino hierarchy. In addition to the declared
fixed point, the package contains a coarse-to-local conditional grid in
`alpha` and `Omega_m0`, with the primordial sector and calibration held fixed.
At its final local grid minimum, a one-dimensional `A_planck` profile is also
evaluated for both the KGB and matched fixed LCDM spectra. This is not an
optimization over all cosmological and nuisance parameters, a posterior,
Bayesian evidence calculation, or a replacement claim for LambdaCDM.

The precise physical status and the empirical boundary are in
`docs/physical_status.md`.
