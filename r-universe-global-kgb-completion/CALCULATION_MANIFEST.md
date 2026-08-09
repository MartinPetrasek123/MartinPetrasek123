# Calculation manifest

## Identity

- Package: `r-universe-global-kgb-completion`
- Theory: globally completed R-alpha luminal kinetic-gravity-braiding action
- Background parameters: `Omega_m0=0.3022734375`, `Omega_r0=9.083909e-5`, `alpha=0.497421875`
- Reproduction command: `bash scripts/run_all.sh`

## Primary scripts

| Script | Result |
| --- | --- |
| `scripts/ru_kgb.py` | Exact global R-alpha background and KGB action functions |
| `scripts/validate_ru_kgb.py` | Background, action, scalar, and tensor gates |
| `scripts/stability_scan.py` | Positivity scan in alpha and Omega_m0 |
| `scripts/matter_qs.py` | Action-derived high-k matter response and growth |
| `scripts/ppn_screening.py` | Cubic-KGB Vainshtein/Cassini gate |
| `scripts/cmb_prerecombination.py` | Pre-recombination background and sound-horizon gate |
| `scripts/make_figures.py` | Stability figure |

## Scope discipline

This manifest records calculations that are actually executed by the command
above. In particular, it does not label a pre-recombination calculation as a
Planck `C_ell` likelihood, a quasistatic matter prediction as a full RSD or
lensing likelihood, or a Vainshtein estimate as an ephemeris fit.
