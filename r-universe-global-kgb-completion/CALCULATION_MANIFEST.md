# Calculation manifest

## Identity

- Package: `r-universe-global-kgb-completion`
- Theory: globally completed R-alpha luminal kinetic-gravity-braiding action
- Background parameters: `Omega_m0=0.3022734375`, `Omega_r0=9.083909e-5`, `alpha=0.497421875`
- Reproduction command: `bash scripts/run_all.sh`
- Native CMB solver: `HEFTCAMB_BIN=/path/to/camb python3 scripts/run_heftcamb_rph.py`
- Fixed Planck point: `CLIPY_SOURCE=/path/to/clipy PLANCK_2018_BASE=/path/to/plc_3.0 /path/to/python-with-clipy scripts/evaluate_planck_2018_fixed.py --spectra-dir generated/heftcamb/convergence/nodes_601`

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

`docs/physical_status.md` records the sharper physical conclusion: this is a
covariant KGB candidate whose fixed-point CMB calculation does not by itself
establish an empirical preference or a derivation from RCD.
