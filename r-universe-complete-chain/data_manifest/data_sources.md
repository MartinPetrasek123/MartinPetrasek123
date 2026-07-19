# Data sources and local mapping

This file records the public data sources used by the R-Universe reproducibility package.

## Pantheon+

- Public repository: https://github.com/PantheonPlusSH0ES/DataRelease
- Paper: Brout et al. 2022, The Pantheon+ Analysis: Cosmological Constraints.
- Local source files used by the fitting code:
  - `pantheon_plus.dat`
  - `pantheon_plus_stat_sys.cov`
- Local historical workspace path:
  - `/Users/mpetr/Documents/Codex/2026-07-09/v/work/r_universe_fit/data/pantheon_plus.dat`
  - `/Users/mpetr/Documents/Codex/2026-07-09/v/work/r_universe_fit/data/pantheon_plus_stat_sys.cov`
- Selection used:
  - remove Cepheid calibrators;
  - require `zHD > 0.01`;
  - resulting sample size: 1580 supernovae.

## DESI BAO

- Public repository used for data vectors and covariance matrices: https://github.com/CobayaSampler/bao_data
- DESI DR1 papers: arXiv:2404.03000, arXiv:2404.03001, arXiv:2404.03002.
- DESI DR2 BAO paper: arXiv:2503.14738.
- Local source files used:
  - `work/public_bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt`
  - `work/public_bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt`
  - historical DESI DR1 files under `/Users/mpetr/Documents/Codex/2026-07-09/v/work/r_universe_fit/data/`.

## Cosmic chronometers

- The current code uses the transparent 31-point compilation embedded in `fit_r_universe.py`.
- Sources represented in the compilation include Simon, Stern, Moresco, Zhang and Ratsimbazafy measurements.
- The values are explicitly listed in the code so that the likelihood can be audited without hidden data files.

## DES-Dovekie / DES-SN5YR

- Public repository: https://github.com/des-science/DES-SN5YR
- DOI: https://doi.org/10.5281/zenodo.12720778
- Local files used:
  - `work/DES-SN5YR/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv`
  - `work/DES-SN5YR/4_DISTANCES_COVMAT/STAT+SYS.npz`
- Use in this package:
  - independent supernova cross-check;
  - not added to Pantheon+ because low-z and calibration/systematics overlap may exist.

## Growth and lensing compressed checks

- RSD base compilation: Sagredo, Nesseris and Sapone 2018.
- Weak-lensing compressed checks:
  - DES Y3 cosmic shear, Amon et al. 2022;
  - KiDS-1000, Asgari et al. 2021;
  - HSC Y3 compressed S8 value used in the local pilot package.

## Important limitation

The BAO fits in this package fix `r_d = 147.09 Mpc`. This package is a complete late-time-background reproducibility package, not yet a full early-universe Boltzmann/CMB implementation.
