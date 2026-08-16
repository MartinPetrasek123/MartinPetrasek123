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
- `generated/kgb_multprobe_conditional_summary.json` - ledger joining the executed fixed-input Planck, Pantheon+, DESI DR2 BAO, chronometer, PPN, and native-RSD-audit records without conflating them into a global posterior;
- `generated/kgb_joint_posterior_production_seeded_v4/` - the released four-chain exact KGB posterior, its production and deep point-cache audits, the immutable execution contract, and the generated manuscript block;
- `generated/rsd_native_best_point/rsd_native_dz_0p005.json` and `generated/rsd_native_best_point_dz_0p0025/rsd_native_dz_0p0025.json` - native H-EFTCAMB `f sigma8` audits and finite-difference convergence records;
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
At the selected matched point `(alpha, Omega_m0)=(0.0975, 0.3075)`, separate
fixed-spectrum `A_planck` profiles and exact-background Pantheon+--DESI DR2
BAO--chronometer evaluations give a conditional KGB-minus-LCDM sum of
`-2.5637653`. This is not a posterior or evidence calculation: the two
`A_planck` values are independently profiled and all remaining
cosmological/nuisance inputs are fixed. The same KGB point has a native linear
`f sigma8` residual audit, but the supplied RSD compilation lacks the survey
covariance, window, AP, and nonlinear ingredients required for a likelihood
and is deliberately excluded from that sum. The KGB action remains a physical
replacement candidate, with LCDM used only as an observational benchmark; no
empirical replacement is asserted.

The precise physical status and the empirical boundary are in
`docs/physical_status.md`.

## Exact posterior release

`generated/kgb_joint_posterior_production_seeded_v4/` is the public,
hash-preserving record for the reported exact four-chain KGB posterior. It
contains the four chains, their saved Cobaya inputs and checkpoints, the
production summary, the independent production-bundle audit, the deep
point-cache audit, and the TeX block used by the manuscript. The release also
contains the exact configuration, proposal covariance, executable sources, and
H-EFTCAMB template recorded by the execution contract.

Verify the released posterior without modifying any artifact:

```sh
PYTHONPATH=. python3 inference/verify_kgb_exact_posterior_release.py
```

The verifier checks the published SHA-256 values, regenerates the posterior
summary from the four released chains, validates the production contract, and
checks the source hashes recorded for the generator, Planck wrapper, late-time
evaluator, and joint evaluator. See
`docs/exact_posterior_release.md` for the precise scope and external
dependencies.
