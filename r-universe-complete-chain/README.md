# R-Universe complete-chain reproducibility package

This package accompanies the manuscript:

**R-Universe as a Testable Extension of Flat LCDM: a Complete Late-Time Expansion, DESI DR2 BAO, DES-Dovekie and Growth-of-Structure Study**

## Contents

- `main.tex` - English manuscript.
- `references.bib` - bibliography.
- `graf.tex` - all figure environments used by the manuscript.
- `figures/` - generated figures.
- `data/raw/` - exact public data files used by the released scripts.
- `tables/all_model_fits.csv` - all fitted model summaries.
- `tables/desi_dr2_bao_predictions.csv` - DESI DR2 BAO predictions and diagnostic pulls.
- `tables/camb_lensed_cls_lcdm_reference.csv` - CAMB reference CMB spectra.
- `tables/camb_matter_power_lcdm_reference.csv` - CAMB reference matter power.
- `code/r_universe_core.py` - shared data loaders, background solver and likelihood blocks.
- `code/extended_fit.py` - Pantheon+ full-covariance + DESI DR1/DR2 + cosmic-chronometer fits.
- `code/des_dovekie_fit.py` - DES-Dovekie STAT+SYS + DESI DR2 + cosmic-chronometer cross-check.
- `code/make_figures_and_tables.py` - figure and table generation.
- `code/derived_predictions.py` - cosmographic, growth and null-test predictions.
- `code/numerical_validation.py` - Newton/bisection, LCDM-limit and distance-integration validation.
- `code/boltzmann_camb.py` - CAMB reference early-universe calculation.
- `code/*.json` - generated numerical results used by the manuscript.
- `data_manifest/data_sources.md` - source URLs and local file mapping.
- `requirements.txt`, `environment.yml` - software environment records.

## Reproduction

Use Python 3.10+ with `numpy`, `scipy`, `pandas` and `matplotlib`.
From this package root:

```bash
python3 code/extended_fit.py
python3 code/des_dovekie_fit.py
python3 code/derived_predictions.py
python3 code/numerical_validation.py
python3 code/boltzmann_camb.py
python3 code/make_figures_and_tables.py
```

The scripts expect the public data files listed in `data_manifest/data_sources.md`.
The generated JSON, CSV and figure files are the numerical source of the manuscript tables and plots.

## Versioning and archive

- Manuscript release tag: `v1.1.0`
- GitHub package: <https://github.com/MartinPetrasek123/MartinPetrasek123/tree/main/r-universe-complete-chain>
- A Zenodo DOI should be minted from the `v1.1.0` GitHub release before journal submission. This draft does not invent or claim a DOI that has not been issued.
- Data filenames, masks, covariance files, vector order and SHA256 checksums are recorded in `data_manifest/data_sources.md`.

## Main numerical result

For Pantheon+ full covariance + DESI DR2 BAO + cosmic chronometers:

- R1 vs LCDM: Delta chi2 = -3.7117, Delta AIC = -1.7117, Delta BIC = +3.6809.

For DES-Dovekie STAT+SYS + DESI DR2 BAO + cosmic chronometers:

- R1 vs LCDM: Delta chi2 = -4.2133, Delta AIC = -2.2133, Delta BIC = +3.3172.

The R1 branch is therefore AIC-favored but not BIC-favored in the present late-time data-only implementation.
