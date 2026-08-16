# R-Universe full manuscript and reproducibility package

This directory is the authoritative source bundle for the manuscript
**Relational Capacity Dynamics, Relational Foliation Gravity, RFG-R, and the
Covariant KGB Realization**.

The release is deliberately explicit about scope. It contains derived
classical RCD/RFG/RFG-R results and a separately specified covariant KGB
realization. The KGB package includes an audited posterior for its stated
Planck 2018 CMB plus lensing, Pantheon+, DESI DR2 BAO, and
cosmic-chronometer probe set. It does not claim a Bayesian-evidence result, a
matched full model comparison with LambdaCDM, an official RSD likelihood, or
an executed CMB/matter likelihood for RFG-R or RFG-R-Xi.

## Authoritative artifacts

- `main.tex` - complete manuscript source.
- `R_Universe_v1.6.1.pdf` - compiled manuscript corresponding to this tag.
- `CITATION.cff` and `codemeta.json` - citation and software metadata.
- `data_manifest/data_sources.md` - input provenance, selections, covariance
  treatment, and SHA-256 fingerprints for the late-time data blocks.
- `rfg-r-completion/r_universe_completion/` - RFG-R and RFG-R-Xi action,
  validation, EFT-map, and constraint-audit source records.
- `../r-universe-global-kgb-completion/` - covariant KGB action,
  H--EFTCAMB contracts, released chains, production audit, and read-only
  posterior verifier.

## Reproduction routes

The historical late-time background scripts and their public input files are
kept in this directory under `code/` and `data/raw/`. Their declared runtime
is Python 3.10 or newer with the packages listed in `requirements.txt`:

```sh
python3 code/run_all.py
```

The KGB exact-posterior release is independently verifiable without a long
solver campaign. From `../r-universe-global-kgb-completion/`, use Python 3.9
or newer with NumPy and SciPy:

```sh
PYTHONPATH=. python3 inference/verify_kgb_exact_posterior_release.py
```

This second command validates the released chain, contract, executable, and
audit hashes, then regenerates the posterior summary from the four released
chains. It does not redistribute the licensed Planck likelihood files or
rerun the 46,206 cached solver evaluations.

## Versioning

- Manuscript release: `v1.6.1`.
- Source snapshot: <https://github.com/MartinPetrasek123/MartinPetrasek123/tree/v1.6.1/r-universe-complete-chain>
- A Zenodo DOI may be minted from this immutable tag before submission. No
  DOI is claimed until one is issued.

The MIT license covers the repository's original code and documentation.
Third-party observational inputs retain their source terms and must be cited
according to `data_manifest/data_sources.md`.
