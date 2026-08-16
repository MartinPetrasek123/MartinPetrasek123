#!/usr/bin/env python3
"""Verify the portable exact-KGB posterior release without altering its records."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUNDLE = ROOT / "generated" / "kgb_joint_posterior_production_seeded_v4"
CHAIN_INDICES = (5, 6, 7, 8)
EXECUTABLE_HASHES = {
    "generator_sha256": ROOT / "scripts" / "generate_heftcamb_rph.py",
    "planck_wrapper_sha256": ROOT / "scripts" / "evaluate_planck_2018_fixed.py",
    "late_time_evaluator_sha256": ROOT / "scripts" / "evaluate_kgb_late_time.py",
    "joint_evaluator_sha256": ROOT / "scripts" / "evaluate_kgb_joint_point.py",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"expected a JSON object: {path}")
    return value


def required_path(bundle: Path, name: str) -> Path:
    path = bundle / name
    if not path.is_file():
        raise FileNotFoundError(f"required release artifact is unavailable: {path}")
    return path


def verify_chain_records(bundle: Path, summary: dict[str, Any], audit: dict[str, Any]) -> None:
    expected_names = [f"chain_{index:02d}.1.txt" for index in CHAIN_INDICES]
    recorded_names = [Path(value).name for value in summary.get("chains", [])]
    if recorded_names != expected_names:
        raise RuntimeError("summary does not identify the four released production chains")
    records = audit.get("chains")
    if not isinstance(records, list) or len(records) != len(CHAIN_INDICES):
        raise RuntimeError("production audit does not contain four chain records")
    for index, record in zip(CHAIN_INDICES, records, strict=True):
        if record.get("index") != index:
            raise RuntimeError("production audit has a mismatched chain index")
        prefix = f"chain_{index:02d}"
        expected = {
            "sample": f"{prefix}.1.txt",
            "input": f"{prefix}.input.yaml",
            "updated": f"{prefix}.updated.yaml",
            "checkpoint": f"{prefix}.checkpoint",
        }
        hashes = record.get("sha256")
        if not isinstance(hashes, dict):
            raise RuntimeError(f"production audit has no hashes for {prefix}")
        for key, name in expected.items():
            if record.get(key) != name or sha256(required_path(bundle, name)) != hashes.get(key):
                raise RuntimeError(f"release hash mismatch for {name}")


def recompute_summary(bundle: Path, stored: dict[str, Any]) -> None:
    chains = [required_path(bundle, f"chain_{index:02d}.1.txt") for index in CHAIN_INDICES]
    with tempfile.TemporaryDirectory(prefix="kgb-posterior-release-") as temporary:
        output = Path(temporary) / "summary.json"
        command = [
            sys.executable,
            str(ROOT / "inference" / "analyze_kgb_joint_chains.py"),
            "--chains",
            *(str(path) for path in chains),
            "--output",
            str(output),
            "--burn-in-fraction",
            str(stored["burn_in_fraction"]),
            "--rhat-threshold",
            str(stored["publication_gate"]["rank_normalized_split_rhat_less_than"]),
            "--minimum-bulk-ess",
            str(stored["publication_gate"]["minimum_bulk_ess"]),
            "--minimum-tail-ess",
            str(stored["publication_gate"]["minimum_tail_ess"]),
        ]
        subprocess.run(command, check=True, capture_output=True, text=True)
        recomputed = load_json(output)
    # The original record intentionally retains its historical absolute paths.
    # Every posterior-derived field must nevertheless be identical after relocation.
    stored_portable = {key: value for key, value in stored.items() if key != "chains"}
    recomputed_portable = {key: value for key, value in recomputed.items() if key != "chains"}
    if recomputed_portable != stored_portable:
        raise RuntimeError("recomputed posterior summary differs from the released record")


def verify(bundle: Path) -> dict[str, Any]:
    summary_path = required_path(bundle, "production_summary.json")
    audit_path = required_path(bundle, "production_audit.json")
    cache_audit_path = required_path(bundle, "point_cache_audit.json")
    summary = load_json(summary_path)
    audit = load_json(audit_path)
    cache_audit = load_json(cache_audit_path)
    if audit.get("status") != "passed" or audit.get("summary_sha256") != sha256(summary_path):
        raise RuntimeError("production audit does not validate the released summary")
    if cache_audit.get("status") != "passed":
        raise RuntimeError("deep point-cache audit is not recorded as passed")
    verify_chain_records(bundle, summary, audit)
    contract_path = required_path(bundle, "campaign_contract.json")
    if audit.get("campaign_contract_sha256") != sha256(contract_path):
        raise RuntimeError("campaign contract hash does not match the production audit")
    threads = required_path(bundle, "threading_verification_omp2.json")
    if audit.get("threading_verification_sha256") != sha256(threads):
        raise RuntimeError("threading verification hash does not match the production audit")
    execution = load_json(required_path(bundle, "execution_contract/manifest.json"))["execution_contract"]
    for field, path in EXECUTABLE_HASHES.items():
        if execution.get(field) != sha256(path):
            raise RuntimeError(f"execution-contract hash mismatch for {path.name}")
    recompute_summary(bundle, summary)
    return {
        "status": "passed",
        "summary_sha256": sha256(summary_path),
        "production_audit_sha256": sha256(audit_path),
        "point_cache_audit_sha256": sha256(cache_audit_path),
        "chains": [f"chain_{index:02d}.1.txt" for index in CHAIN_INDICES],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    args = parser.parse_args()
    print(json.dumps(verify(args.bundle.resolve()), indent=2))


if __name__ == "__main__":
    main()
