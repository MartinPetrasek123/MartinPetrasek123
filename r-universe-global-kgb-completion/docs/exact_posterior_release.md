# Exact KGB Posterior Release

## Contents

This versioned release makes the reported KGB posterior independently
verifiable without a long new inference run. The release includes:

- the four original Cobaya chain files;
- their immutable input files, updated inputs, checkpoints, and sampler
  covariances;
- `production_summary.json`, `production_audit.json`, and
  `point_cache_audit.json`;
- the campaign contract, threading check, execution-contract manifest,
  production configuration, initial proposal covariance, and spectrum
  template;
- the posterior rendering source `posterior_results.tex`;
- the runner, analyzer, production-bundle auditor, point-cache auditor, joint
  likelihood wrapper, and exact source files whose hashes were recorded in the
  execution contract.

The numerical result is the posterior for the stated covariant KGB action and
the stated Planck 2018 CMB plus lensing, Pantheon+, DESI DR2 BAO, and
cosmic-chronometer probe set. It is not a Bayesian evidence calculation, an
RSD likelihood analysis, or a model comparison with LambdaCDM.

## Verification

From the package root, run:

```sh
PYTHONPATH=. python3 inference/verify_kgb_exact_posterior_release.py
```

The portable verifier requires Python 3.9 or newer and NumPy plus SciPy. It
accepts only platform-level floating-point rounding in regenerated summary
values (relative or absolute tolerance \(10^{-12}\)); all released files,
contract records, and executable inputs are verified by exact SHA-256 hashes.
Its scope is a deterministic validation of released artifacts, not a fresh
cosmological sampling campaign.

The command is read-only. It validates the production summary and contract
hashes, the four chain/input/updated/checkpoint hash records, the threading
record, the executable hashes in the execution-contract manifest, and the
deep point-cache audit status. It then regenerates all posterior-derived
summary fields from the released chains and requires them to agree with the
published summary.

The historical summary records absolute paths from its execution host. The
release verifier deliberately compares their filenames after relocation while
requiring every content hash to match. It never rewrites the historical
summary, audit, or chain files.

## External Dependencies

Re-executing the full likelihood, rather than verifying the released posterior,
requires Cobaya, NumPy, SciPy, PyYAML, a matching H-EFTCAMB build, the
official Planck 2018 likelihood distribution under its own license, and the
observational inputs specified by the production configuration. The release
preserves the relevant source and binary/template hashes, but does not
redistribute licensed Planck data, the external solver binary, or the
46,206 per-point cache archives. The latter were independently checked before
the recorded deep point-cache audit was issued.
