# Executable KGB Posterior

`kgb_joint_planck_late.yaml` specifies an eight-parameter Cobaya chain.  Each
likelihood call invokes `scripts/evaluate_kgb_joint_point.py`, which regenerates
the exact R-Universe covariant-KGB RPH functions, runs H-EFTCAMB, then applies
the official Planck 2018 Plik-lite plus low-ell plus lensing likelihoods and
the exact-background Pantheon+, DESI DR2 BAO, and chronometer likelihoods.
The Planck contribution explicitly includes its documented absolute-calibration
prior, \(A_{\rm Planck}=1\pm0.0025\); it is not silently supplied by the basic
CLIPY likelihood object.
Completed per-point solver artifacts use a verified archive manifest. The
`full` profile is lossless; production uses the explicit `likelihood-core`
profile, which retains every likelihood-defining input, spectrum, report, log,
and executable hash and fingerprints each omitted recomputable diagnostic
output. Neither profile changes the likelihood ordinate.
On a host where the configured interpreter is cloud-backed, set
`KGB_EVALUATOR_PYTHON` to an identical local interpreter copy. This only
changes process placement; the declared action, parameters, data, and
likelihood code remain the same.

Run from this package root with the declared interpreter:

```sh
PYTHONPATH=. /path/to/python inference/run_kgb_joint_chain.py \
  --output generated/kgb_joint_chain/chain_01 --max-samples 4000 --seed 101
```

The immutable final-production target is
`kgb_joint_planck_late_production_final.yaml`. Its isolated cache key includes
SHA-256 fingerprints of the joint evaluator, RPH generator, Planck wrapper,
late-time evaluator, H-EFTCAMB binary, and spectrum template. Before any
manuscript update, `snapshot_kgb_execution_contract.py` preserves one
digest-verified campaign-level copy of the configuration and shared spectrum
template; each point keeps its own independently verified archive. The
parallel controller records this contract, starts four independently seeded
chains in two-solver waves, and updates the desktop manuscript only after the
external convergence and ESS gate passes:

```sh
PYTHONPATH=. /path/to/python inference/run_kgb_joint_parallel_campaign.py \
  --config inference/kgb_joint_planck_late_production_final.yaml \
  --output-directory generated/kgb_joint_posterior_production_final \
  --initial-states 1000 --extension-states 1000 --maximum-states 4000 \
  --maximum-parallel 2 --seed-base 290000
```

For an already started two-slot campaign on a host with spare physical cores,
`accelerate_kgb_joint_campaign.py` can add the missing chains 07 and 08 while
the inherited 05 and 06 processes finish. It writes `accelerated_execution.json`
and preserves the exact configuration, likelihood contract, numerical nodes,
and independent seeds; only process scheduling changes. It must never merge an
earlier cache that lacks the current execution provenance.

## Production posterior protocol

One chain is not a posterior result. Four independently seeded chains run in
two independent solver slots, a limit chosen to avoid oversubscribing the
native solver. A locally estimated Hessian covariance is only an MCMC proposal; it
does not supply an uncertainty estimate.  The supplied production helpers
refuse to overwrite an existing chain and reject chains that stop below the
declared retained-weight threshold:

```sh
PYTHONPATH=. /path/to/python inference/run_kgb_joint_production_ensemble.py \
  --max-samples 1000 --seed-base 600
```

After the four files `chain_05.1.txt` through `chain_08.1.txt` exist, summarize
them without modifying samples:

```sh
PYTHONPATH=. /path/to/python inference/analyze_kgb_joint_chains.py \
  --chains generated/kgb_joint_posterior_calibrated/chain_05.1.txt \
           generated/kgb_joint_posterior_calibrated/chain_06.1.txt \
           generated/kgb_joint_posterior_calibrated/chain_07.1.txt \
           generated/kgb_joint_posterior_calibrated/chain_08.1.txt \
  --burn-in-fraction 0.25 \
  --output generated/kgb_joint_posterior_calibrated/production_summary.json
```

Only a summary whose rank-normalized, folded split-\(\hat R\) is below 1.01 for
every sampled parameter and whose rank-normalized bulk and 5\%/95\% tail ESS
are each at least 400 may be described as converged. The simple split-\(\hat R\)
is retained in the JSON output as an audit diagnostic. The minimum likelihood
retained by the chain is not a Bayesian evidence and does not establish a
comparison with another cosmology. The evaluated probe set is Planck 2018 CMB plus lensing,
Pantheon+, DESI DR2 BAO, and cosmic chronometers.

Before a converged summary is rendered into the manuscript,
`audit_kgb_production_bundle.py` independently checks the four saved chain
contracts, recomputes their retained weights, checks the recorded publication
gate, and records SHA-256 hashes for the chain, input, updated, and checkpoint
files. Proposal scale and stopping count may vary between independent chains;
the likelihood definition, priors, and target sampler controls may not.
`finalize_kgb_production.py` reruns those guards, including the completed-point
cache audit and the matching shared execution-contract snapshot, before it
writes the final TeX block into the desktop manuscript.
`wait_for_kgb_production_finalization.py` can supervise that finalization for a
long-running controller without treating an unconverged intermediate summary
as publishable.

`wait_and_run_kgb_joint_ensemble.py` is the bounded production supervisor used
when an initial chain is already running. It starts the remaining three
chains, writes the common summary, and, only if the split-\(\hat R\) gate fails,
resumes all four compatible Cobaya checkpoints to a declared larger stored-
state limit. Rejected proposals remain represented by the standard Cobaya
sample weights. A resumed sampler retains its saved proposal state; the target
posterior, data, and priors are unchanged. The runner reloads that chain's
saved ``.input.yaml`` when resuming, so later edits to the template cannot
silently alter the historical target; only the declared stopping limit may be
extended.

The RSD compilation is deliberately excluded.  It does not contain the
survey-complete covariance, window, Alcock-Paczynski mapping, and nonlinear
nuisance likelihood definitions required for a valid RSD likelihood.

## Point-cache archival

Each exact likelihood point can contain several solver spectra and logs.  To
avoid allowing the immutable audit trail to exhaust the local volume,
`archive_kgb_joint_cache.py` creates an xz-compressed tar archive for every completed
point, verifies each retained member's name, size, and SHA-256 content hash
before pruning, records SHA-256 hashes for the archive and every retained
member, and leaves `joint_point.json` in place for exact cache reuse. The
`full` policy preserves every raw artifact. The `likelihood-core` policy also
records the SHA-256 and size of the intentionally omitted, recomputable
diagnostic files, while preserving all likelihood-defining artifacts.
Only the explicit `--prune` mode removes raw files, and only after this
verification succeeds. Its batch command skips points with a completed summary
newer than five minutes by default, avoiding any race with a running evaluator:

```sh
PYTHONPATH=. /path/to/python inference/archive_kgb_joint_cache.py \
  --cache-directory generated/kgb_joint_chain_points_calibrated --prune
```

No likelihood datum is discarded. Restore and hash-check retained artifacts
before inspection with `inference/restore_kgb_joint_cache.py --point-directory
<point-directory>`.
