# Calculation manifest

## Identity

- Package: `r-universe-global-kgb-completion`
- Theory: globally completed R-alpha luminal kinetic-gravity-braiding action
- Background parameters: `Omega_m0=0.3022734375`, `Omega_r0=9.083909e-5`, `alpha=0.497421875`
- Reproduction command: `bash scripts/run_all.sh`
- Native CMB solver: `HEFTCAMB_BIN=/path/to/camb python3 scripts/run_heftcamb_rph.py`
- Fixed Planck point: `CLIPY_SOURCE=/path/to/clipy PLANCK_2018_BASE=/path/to/plc_3.0 /path/to/python-with-clipy scripts/evaluate_planck_2018_fixed.py --spectra-dir generated/heftcamb/convergence/nodes_601`
- Conditional KGB profile: `scripts/profile_planck_kgb.py --binary /path/to/camb --template /path/to/HighLExtrapTemplate_lenspotentialCls.dat --clipy-source /path/to/clipy --planck-base /path/to/plc_3.0 --python /path/to/python-with-clipy --alpha-values 0.030 0.035 0.040 0.045 --omega-m0-values 0.30825 0.30850 0.30875 --output generated/planck_profile_final_local`
- Fixed-spectrum calibration profiles: `scripts/profile_planck_calibration.py --spectra-dir /path/to/spectra --clipy-source /path/to/clipy --planck-base /path/to/plc_3.0 --python /path/to/python-with-clipy --a-planck-values 1.0010 1.0015 1.0020 1.0025 1.0030 --output /path/to/profile-output`
- Exact late-time geometry: `python3 scripts/evaluate_kgb_late_time.py --data-root /path/to/r-universe-complete-chain/data/raw --solver-dir /path/to/executed-solver-point --omega-m0 0.3075 --omega-r0 9.083909e-5 --alpha 0.0975 --integration-nodes 32769 --output /path/to/result.json`
- Native RSD audit: generate multi-redshift spectra with `--transfer-redshifts`, then run `python3 scripts/evaluate_kgb_rsd_native.py --solver-dir /path/to/executed-solver-point --rsd-data /path/to/rsd.csv --finite-difference-dz 0.005 --output /path/to/result.json`

## Primary scripts

| Script | Result |
| --- | --- |
| `scripts/ru_kgb.py` | Exact global R-alpha background and KGB action functions |
| `scripts/validate_ru_kgb.py` | Background, action, scalar, and tensor gates |
| `scripts/validate_independent_kgb_identities.py` | Independent derivative, KGB alpha-function, and scalar-equation identities |
| `scripts/stability_scan.py` | Positivity scan in alpha and Omega_m0 |
| `scripts/matter_qs.py` | Action-derived high-k matter response and growth |
| `scripts/ppn_screening.py` | Cubic-KGB Vainshtein/Cassini gate |
| `scripts/cmb_prerecombination.py` | Pre-recombination background and sound-horizon gate |
| `scripts/generate_heftcamb_rph.py` | Exact KGB functions in the native H-EFTCAMB RPH convention |
| `scripts/run_heftcamb_rph.py` | Native CMB and linear-matter solver plus spline convergence |
| `scripts/evaluate_planck_2018_fixed.py` | Official Planck 2018 fixed-point CMB and lensing likelihood |
| `scripts/summarize_planck_fixed_runs.py` | Fixed-point likelihood and numerical-convergence audit record |
| `scripts/profile_planck_kgb.py` | Fixed-input local KGB Planck grid with a matched fixed LCDM reference |
| `scripts/profile_planck_calibration.py` | One-dimensional Planck absolute-calibration grid at fixed spectra |
| `scripts/evaluate_kgb_late_time.py` | Pantheon+ full-covariance, DESI DR2 BAO full-covariance, and cosmic-chronometer geometric likelihood |
| `scripts/assemble_kgb_joint_profile.py` | Adds executed fixed-input Planck ordinates to the exact late-time likelihood |
| `scripts/merge_kgb_joint_profiles.py` | Deduplicates conditional joint-grid chunks without altering point values |
| `scripts/evaluate_kgb_rsd_native.py` | Direct `sigma8` integration and native `f sigma8` residual audit from H-EFTCAMB `P(k,z)` |
| `scripts/summarize_kgb_multprobe.py` | Machine-readable ledger of the executed conditional records |
| `scripts/make_figures.py` | Stability figure |

## Scope discipline

`bash scripts/run_all.sh` executes the action-level gates. The native
H-EFTCAMB and Planck steps require the external solver and the official
likelihood distribution, so their already executed reference outputs are
stored separately under `generated/heftcamb/` and
`generated/planck_2018_fixed_summary.json`. The latter is an actual Planck
`C_ell` and lensing likelihood evaluation at one fixed point. It is not a
full RSD/BAO likelihood, ephemeris fit, parameter posterior, evidence
calculation, or model comparison.

`generated/kgb_multprobe_conditional_summary.json` records a selected matched
KGB--LCDM pair at `alpha=0.0975`, `Omega_m0=0.3075`: its own executed
fixed-spectrum `A_planck` minimum is added to an exact-background Pantheon+
full-covariance, DESI DR2 BAO full-covariance, and chronometer statistic. The
reported KGB-minus-LCDM conditional sum is `-2.563765256550596`. This is not a
posterior, evidence calculation, or an optimization over cosmological and
nuisance parameters. The included native RSD output is an intentionally
separate residual audit because the local compilation lacks the full survey
covariance/window/AP/nonlinear likelihood inputs.

`docs/physical_status.md` records the sharper physical conclusion: this is a
covariant KGB candidate whose fixed-point and limited conditional CMB
calculations do not establish an empirical preference or a derivation from
RCD.
