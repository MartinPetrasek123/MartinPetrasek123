#!/usr/bin/env python3
"""Independently audit a converged four-chain KGB production bundle.

This is a publication guard, not a likelihood calculation. It verifies that
all chains share the same physical target, recomputes retained weights, checks
the stored convergence gate, and records file hashes for the posterior inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import yaml

from inference.analyze_kgb_joint_chains import read_chain


ROOT = Path(__file__).resolve().parents[1]
CHAIN_INDICES = (5, 6, 7, 8)
TARGET_SAMPLER_KEYS = (
    "covmat",
    "learn_proposal",
    "learn_every",
    "Rminus1_stop",
    "Rminus1_cl_stop",
    "max_tries",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_target(input_path: Path) -> dict[str, Any]:
    """Extract posterior-defining settings and independent chain identities."""
    loaded = yaml.safe_load(input_path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict):
        raise RuntimeError(f"chain input is not a mapping: {input_path}")
    try:
        sampler = loaded["sampler"]["mcmc"]
        return {
            "likelihood": loaded["likelihood"],
            "params": loaded["params"],
            "sampler_target_controls": {key: sampler[key] for key in TARGET_SAMPLER_KEYS},
            "sampler_seed": sampler.get("seed"),
        }
    except (KeyError, TypeError) as error:
        raise RuntimeError(f"chain input is incomplete: {input_path}") from error


def production_gate(summary: dict[str, Any]) -> tuple[float, float, float]:
    gate = summary.get("publication_gate")
    if not isinstance(gate, dict):
        raise RuntimeError("summary lacks a recorded publication gate")
    rhat = float(gate["rank_normalized_split_rhat_less_than"])
    bulk = float(gate["minimum_bulk_ess"])
    tail = float(gate["minimum_tail_ess"])
    if (
        not 1.0 < rhat <= 1.1
        or bulk <= 0.0
        or tail <= 0.0
        or not all(math.isfinite(value) for value in (rhat, bulk, tail))
    ):
        raise RuntimeError("summary has an invalid publication gate")
    return rhat, bulk, tail


def assert_gate(summary: dict[str, Any], rhat: float, bulk: float, tail: float) -> None:
    if summary.get("passes_publication_gate") is not True:
        raise RuntimeError("summary does not pass its recorded publication gate")
    parameters = summary.get("parameters")
    if not isinstance(parameters, dict) or not parameters:
        raise RuntimeError("summary lacks parameter diagnostics")
    for name, diagnostic in parameters.items():
        if not isinstance(diagnostic, dict):
            raise RuntimeError(f"invalid diagnostic record: {name}")
        try:
            values = (
                float(diagnostic["rank_normalized_split_rhat"]),
                float(diagnostic["bulk_ess"]),
                float(diagnostic["tail_ess"]),
            )
        except (KeyError, TypeError, ValueError) as error:
            raise RuntimeError(f"incomplete diagnostic record: {name}") from error
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError(f"non-finite diagnostic record: {name}")
        if values[0] >= rhat or values[1] < bulk or values[2] < tail:
            raise RuntimeError(f"publication gate failed for {name}")


def campaign_contract(directory: Path) -> dict[str, Any] | None:
    """Load and validate the additional immutable requirements of a campaign."""
    path = directory / "campaign_contract.json"
    if not path.is_file():
        return None
    loaded = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict):
        raise RuntimeError("campaign contract is not a mapping")
    return loaded


def assert_campaign_execution_record(
    contract: dict[str, Any], directory: Path, seeds: list[Any]
) -> dict[str, Any]:
    """Verify native chain identities and the exact-solver threading record."""
    try:
        expected = contract["sampler_rng"]["chain_seeds"]
        configured_threads = int(contract["omp_num_threads_per_solver"])
        verification = contract["threading_verification"]
        verification_sha256 = str(verification["sha256"])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError("campaign contract lacks required seed or threading provenance") from error
    expected_seeds = [int(expected[str(index)]) for index in CHAIN_INDICES]
    try:
        actual_seeds = [int(seed) for seed in seeds]
    except (TypeError, ValueError) as error:
        raise RuntimeError("production campaign has no native sampler seed for every chain") from error
    if actual_seeds != expected_seeds or len(set(actual_seeds)) != len(actual_seeds):
        raise RuntimeError("production chain seeds do not match the immutable campaign contract")

    threading_path = directory / f"threading_verification_omp{configured_threads}.json"
    if not threading_path.is_file() or sha256(threading_path) != verification_sha256:
        raise RuntimeError("threading-verification artifact is missing or differs from the campaign contract")
    threading = json.loads(threading_path.read_text(encoding="utf-8"))
    try:
        difference = float(threading["benchmark"]["absolute_difference"])
        tolerance = float(threading["benchmark"]["tolerance"])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError("threading-verification artifact is malformed") from error
    if (
        threading.get("passes") is not True
        or threading.get("config_sha256") != contract.get("config_sha256")
        or threading.get("reference_omp_num_threads") != 1
        or threading.get("selected_omp_num_threads") != configured_threads
        or threading.get("openblas_num_threads") != 1
        or not math.isfinite(difference)
        or not math.isfinite(tolerance)
        or difference > tolerance
    ):
        raise RuntimeError("threading-verification artifact does not validate the execution contract")
    return {
        "campaign_contract": "campaign_contract.json",
        "campaign_contract_sha256": sha256(directory / "campaign_contract.json"),
        "threading_verification": threading_path.name,
        "threading_verification_sha256": verification_sha256,
        "omp_num_threads_per_solver": configured_threads,
    }


def audit(summary_path: Path, output_path: Path) -> dict[str, Any]:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if not isinstance(summary, dict):
        raise RuntimeError("posterior summary is not a mapping")
    rhat, bulk, tail = production_gate(summary)
    assert_gate(summary, rhat, bulk, tail)

    burn_fraction = float(summary.get("burn_in_fraction"))
    chains = summary.get("chains")
    if not isinstance(chains, list) or len(chains) != len(CHAIN_INDICES):
        raise RuntimeError("summary does not name exactly four production chains")
    expected_paths = [summary_path.parent / f"chain_{index:02d}.1.txt" for index in CHAIN_INDICES]
    if [Path(path).resolve() for path in chains] != [path.resolve() for path in expected_paths]:
        raise RuntimeError("summary chain list does not match the production bundle")

    targets: list[dict[str, Any]] = []
    chain_records: list[dict[str, Any]] = []
    retained_weights: list[int] = []
    for index, chain_path in zip(CHAIN_INDICES, expected_paths, strict=True):
        prefix = chain_path.parent / chain_path.name.removesuffix(".1.txt")
        input_path = prefix.with_suffix(".input.yaml")
        updated_path = prefix.with_suffix(".updated.yaml")
        checkpoint_path = prefix.with_suffix(".checkpoint")
        for path in (chain_path, input_path, updated_path, checkpoint_path):
            if not path.is_file():
                raise FileNotFoundError(f"production artifact unavailable: {path}")
        _, total_weight, retained_weight = read_chain(chain_path, burn_fraction)
        targets.append(canonical_target(input_path))
        retained_weights.append(retained_weight)
        chain_records.append(
            {
                "index": index,
                "sample": chain_path.name,
                "input": input_path.name,
                "updated": updated_path.name,
                "checkpoint": checkpoint_path.name,
                "total_weight": total_weight,
                "retained_weight": retained_weight,
                "sha256": {
                    "sample": sha256(chain_path),
                    "input": sha256(input_path),
                    "updated": sha256(updated_path),
                    "checkpoint": sha256(checkpoint_path),
                },
            }
        )
    posterior_targets = [
        {key: value for key, value in target.items() if key != "sampler_seed"}
        for target in targets
    ]
    if any(target != posterior_targets[0] for target in posterior_targets[1:]):
        raise RuntimeError("production chains do not share one posterior target")
    seeds = [target["sampler_seed"] for target in targets]
    if any(seed is None for seed in seeds):
        raise RuntimeError("production chains do not record a native sampler seed for every chain")
    if len({int(seed) for seed in seeds}) != len(seeds):
        raise RuntimeError("production chains do not have independent native sampler seeds")
    if summary.get("per_chain_total_weight") != [record["total_weight"] for record in chain_records]:
        raise RuntimeError("summary total weights do not match chain files")
    if summary.get("per_chain_retained_weight") != retained_weights:
        raise RuntimeError("summary retained weights do not match chain files")
    if int(summary.get("combined_retained_weight")) != sum(retained_weights):
        raise RuntimeError("summary combined retained weight does not match chain files")

    contract = campaign_contract(summary_path.parent)
    execution_record = (
        assert_campaign_execution_record(contract, summary_path.parent, seeds)
        if contract is not None
        else {}
    )
    record = {
        "status": "passed",
        "model": summary.get("model"),
        "summary": summary_path.name,
        "summary_sha256": sha256(summary_path),
        "publication_gate": summary["publication_gate"],
        "burn_in_fraction": burn_fraction,
        "target_sampler_controls": posterior_targets[0]["sampler_target_controls"],
        "chain_sampler_seeds": seeds,
        "allowed_chain_specific_sampler_settings": ["output", "proposal_scale", "max_samples", "seed", "checkpoint state"],
        "chains": chain_records,
        **execution_record,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary",
        type=Path,
        default=ROOT / "generated" / "kgb_joint_posterior_calibrated" / "production_summary.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "generated" / "kgb_joint_posterior_calibrated" / "production_audit.json",
    )
    args = parser.parse_args()
    print(json.dumps(audit(args.summary.resolve(), args.output.resolve()), indent=2))


if __name__ == "__main__":
    main()
