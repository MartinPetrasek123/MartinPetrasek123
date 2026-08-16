# Data sources, checksums and calculation provenance

This manifest is part of the reproducibility record for the R-Universe
complete-chain package. It records the public data sources, exact files,
selection masks, covariance usage and code paths used to generate every
number in the manuscript.

## Package release

- Public package: https://github.com/MartinPetrasek123/MartinPetrasek123/tree/main/r-universe-complete-chain
- Manuscript release tag: `v1.6.1`
- Permanent archive: Zenodo DOI must be minted from the `v1.6.1` GitHub
  release before journal submission. No DOI is invented in the draft.

## Pantheon+

- Full dataset name: Pantheon+SH0ES DataRelease, Pantheon+ distance data.
- Authors: Pantheon+SH0ES Collaboration / Brout et al.
- DOI: not explicitly provided in the GitHub data repository; cite Brout et al. 2022 and the public repository.
- Public repository: https://github.com/PantheonPlusSH0ES/DataRelease
- Release/version: public repository main branch; source repository last public update visible as 2022-12-21.
- Download/access date: 2026-07-21.
- License: not explicitly stated in the source repository at access time.
- Reference paper: Brout et al. 2022, The Pantheon+ Analysis: Cosmological
  Constraints.
- Analysis script: `code/extended_fit.py`
- Input table used by the fit: `data/raw/pantheon_plus/pantheon_plus.dat`
- Covariance used by the fit: `data/raw/pantheon_plus/pantheon_plus_stat_sys.cov`
- SHA256:
  - `pantheon_plus.dat`: `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8`
  - `pantheon_plus_stat_sys.cov`: `abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc`
- Selection:
  - remove Cepheid calibrators;
  - require `zHD > 0.01`;
  - retained sample size: 1580 supernovae.
- Calculation path:
  1. read data table and covariance;
  2. apply the same mask to data vector and covariance rows/columns;
  3. compute `D_C(zHD)` and `D_L = (1 + zHEL) D_C(zHD)`;
  4. compute theoretical distance modulus;
  5. analytically marginalize the additive SN intercept;
  6. return the full-covariance `chi2_SN`.

## DESI DR2 BAO

- Full dataset name: DESI DR2 Gaussian BAO all-tracer combined vector and covariance.
- Authors: DESI Collaboration; Cobaya BAO data repository maintainers.
- DOI: DESI DR2 BAO paper DOI is cited through the publication record; repository DOI not explicitly stated.
- Public repository for vectors/covariances: https://github.com/CobayaSampler/bao_data
- Release/version: `bao_data` release `v2.6`, latest public release listed 2025-03-20.
- Download/access date: 2026-07-21.
- License: not explicitly stated in the source repository at access time.
- DESI DR2 BAO paper: arXiv:2503.14738.
- Analysis script: `code/extended_fit.py`
- Input mean vector:
  `data/raw/desi_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt`
- Input covariance:
  `data/raw/desi_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt`
- SHA256:
  - DR2 mean: `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`
  - DR2 covariance: `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`
- Data-vector order used in `tables/desi_dr2_bao_predictions.csv`:
  - `DV_over_rs` at `z = 0.295`;
  - `DM_over_rs`, `DH_over_rs` at `z = 0.510`;
  - `DM_over_rs`, `DH_over_rs` at `z = 0.706`;
  - `DM_over_rs`, `DH_over_rs` at `z = 0.934`;
  - `DM_over_rs`, `DH_over_rs` at `z = 1.321`;
  - `DM_over_rs`, `DH_over_rs` at `z = 1.484`;
  - `DH_over_rs`, `DM_over_rs` at `z = 2.330`.
- Calculation path:
  1. read mean vector and full covariance;
  2. compute `D_M`, `D_H`, `D_V` from the same `E(z)` as the SN fit;
  3. fix `r_d = 147.09 Mpc`;
  4. evaluate the prediction vector in the published order;
  5. compute the full-covariance BAO quadratic form.

## DESI DR1 BAO

- Full dataset name: DESI DR1 Gaussian BAO public likelihood vector and covariance.
- Authors: DESI Collaboration; Cobaya BAO data repository maintainers.
- DOI: DESI DR1 BAO papers cited through publication records; repository DOI not explicitly stated.
- Public repository for vectors/covariances: https://github.com/CobayaSampler/bao_data
- Release/version: `bao_data` public release series, DESI 2024 files.
- Download/access date: 2026-07-21.
- License: not explicitly stated in the source repository at access time.
- DESI DR1 papers: arXiv:2404.03000, arXiv:2404.03001, arXiv:2404.03002.
- Analysis script: `code/extended_fit.py`
- Input mean vector: `data/raw/desi_dr1/desi_2024_bao_mean.txt`
- Input covariance: `data/raw/desi_dr1/desi_2024_bao_cov.txt`
- SHA256:
  - DR1 mean: `dd2873a0b88459a491af3c0c0307ba059f62df9211d5b976760f310565a1be68`
  - DR1 covariance: `bbafa9074b51cf1a45e0d10e4f37db8c0e80a5d1d1788857abb7fc49fb21abcc`
- Use in manuscript:
  - robustness comparison against the main DESI DR2 BAO fit.

## DES-Dovekie / DES-SN5YR

- Full dataset name: Dark Energy Survey Supernova 5YR Data Release / DES-Dovekie distances.
- Authors: DES Collaboration, DES Supernova Program, Sanchez et al., Vincenzi et al.
- Public repository: https://github.com/des-science/DES-SN5YR
- Zenodo DOI: https://doi.org/10.5281/zenodo.12720778
- Release/version: source README notes the original DES-SN5YR results at repository tag `1.3`; this package uses the public files accessed on 2026-07-21.
- Download/access date: 2026-07-21.
- License: not explicitly stated in the source repository at access time.
- Reference papers: DES-SN5YR cosmology release and Sanchez et al. 2024.
- Analysis script: `code/des_dovekie_fit.py`
- Input table: `data/raw/des_dovekie/DES-Dovekie_HD.csv`
- Covariance: `data/raw/des_dovekie/STAT+SYS.npz`
- Metadata audit file: `data/raw/des_dovekie/DES-Dovekie_Metadata.csv`
- SHA256:
  - `DES-Dovekie_HD.csv`: `2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf`
  - `STAT+SYS.npz`: `ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b`
  - `DES-Dovekie_Metadata.csv`: `45ad71f8470eaecfe2b386699ef66b26b0717c50f445d5f32941988d32c75388`
- Use in manuscript:
  - independent SN cross-check;
  - not combined with Pantheon+ because low-redshift anchors and calibration
    systematics may overlap.
- Calculation path:
  1. read DES-Dovekie Hubble diagram and STAT+SYS covariance;
  2. preserve published ordering;
  3. compute distance moduli from the candidate background;
  4. analytically marginalize the additive SN intercept;
  5. combine with DESI DR2 BAO and chronometers.

## Cosmic chronometers

- Full dataset name: 31-point cosmic chronometer H(z) compilation used for transparent late-time stress testing.
- Authors: compiled from Simon, Stern, Moresco, Zhang, Ratsimbazafy and related cosmic-chronometer measurements.
- Reference record: the row-level `reference` column identifies the original
  measurements. The manuscript cites the 31-point compilation of
  Vagnozzi, Loeb and Moresco (arXiv:2011.11645); original-source provenance
  is retained in the CSV rather than being represented as a survey likelihood.
- URL: values are stored in this release at `data/raw/chronometers/cosmic_chronometers.csv`.
- Release/version: package release `v1.6.1`.
- Download/access date: inherited from the previous local pilot compilation; package creation date 2026-07-21.
- License: not applicable as a manually assembled literature-value table; cite original measurements.
- Analysis scripts:
  - `code/extended_fit.py`
  - `code/des_dovekie_fit.py`
- Data representation:
  - the 31 rows are stored in `data/raw/chronometers/cosmic_chronometers.csv`.
- SHA256:
  - `cosmic_chronometers.csv`: `784f6d0387acadb5335f607a2182f48b138eebfa11aa844ce8f9e1904076b36a`
- Sources represented:
  - Simon, Stern, Moresco, Zhang, Ratsimbazafy and later chronometer
    measurements.
- Covariance treatment:
  - diagonal likelihood;
  - no correlated chronometer systematic covariance is assumed in this release.
- Calculation path:
  1. evaluate `H_model(z_i) = H0 E(z_i)`;
  2. compute the diagonal sum over `(H_i - H_model)^2 / sigma_i^2`;
  3. add to SN and BAO blocks.

## Growth and lensing compressed checks

- Diagnostic script: `code/derived_predictions.py`
- RSD audit file:
  `data/raw/growth/rsd_fsigma8_sagredo2018_plus_boss2024.csv`
- Weak-lensing compressed audit file:
  `data/raw/growth/lensing_s8_compressed.csv`
- SHA256:
  - RSD file: `fe6b4e2e488991751e24d21d3c0a9065a5b1767b1d56cf3fae01c157f95f6e4d`
  - lensing file: `f2258e1c7b87241190bbb5bb32d1ea01c4ed992d30374820357a541df33e14f2`
- Use in manuscript:
  - diagnostic consistency only;
  - not claimed as a full RSD/weak-lensing likelihood.

## Output provenance

| Output | Generated by | Numerical inputs |
|---|---|---|
| `code/extended_results.json` | `code/extended_fit.py` | Pantheon+, DESI DR1/DR2 BAO, chronometers |
| `code/des_dovekie_results.json` | `code/des_dovekie_fit.py` | DES-Dovekie, DESI DR2 BAO, chronometers |
| `code/derived_predictions.json` | `code/derived_predictions.py` | best-fit vectors from `extended_results.json` |
| `code/numerical_validation.json` | `code/numerical_validation.py` | independent bisection, Newton solver, LCDM analytic branch |
| `code/boltzmann_camb_results.json` | `code/boltzmann_camb.py` | CAMB LCDM reference early-universe sector |
| `code/robustness_results.json` | `code/robustness_suite.py` | SN-only, BAO-only, CC-only, no-CC, free-rd and BAO jackknife fits |
| `code/profile_likelihood_results.json` | `code/profile_likelihood.py` | theta profile likelihood with H0/Omega_m re-optimization |
| `tables/all_model_fits.csv` | `code/make_figures_and_tables.py` | JSON result files |
| `tables/robustness_suite.csv` | `code/robustness_suite.py` | dataset/sound-horizon robustness summary |
| `tables/profile_likelihood_theta.csv` | `code/profile_likelihood.py` | profiled theta grid |
| `tables/desi_dr2_bao_predictions.csv` | `code/make_figures_and_tables.py` | DESI DR2 vector/covariance and best fits |
| `tables/camb_lensed_cls_lcdm_reference.csv` | `code/boltzmann_camb.py` | CAMB TT/EE/BB/TE lensed spectra |
| `tables/camb_matter_power_lcdm_reference.csv` | `code/boltzmann_camb.py` | CAMB matter power at z=0,0.5,1,2 |
| `figures/*.pdf`, `figures/*.png` | `code/make_figures_and_tables.py` | JSON/CSV result files and best-fit functions |

## Important limitations

- Late-time BAO fits fix `r_d = 147.09 Mpc`; `code/boltzmann_camb.py`
  now computes an LCDM reference value `r_drag = 146.815 Mpc`.
- This release is a complete late-time background and minimal-growth
  reproducibility package.
- It is not yet a full R1 early-universe Boltzmann/CMB likelihood.
- A final publication release should be archived on Zenodo with a DOI.
