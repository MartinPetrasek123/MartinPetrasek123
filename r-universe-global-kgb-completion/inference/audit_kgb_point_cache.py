#!/usr/bin/env python3
"""Audit every completed exact-point record in a KGB likelihood cache.

This verifies the immutable identity of an exact evaluator point independently
of Cobaya: the cache-directory key, component-code hashes, finite component
sum, and verified artifact archive must all match the declared likelihood
contract.  In-progress directories are deliberately ignored, so the command
is also useful as a live integrity check during a production campaign.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from inference.archive_kgb_joint_cache import ARCHIVE, MANIFEST, sha256, verify_archive


LIKELIHOOD_NAME = "likelihoods.kgb_joint_likelihood.KGBJointLikelihood"
PROVENANCE_COMPONENTS = {
    "generator_sha256": "scripts/generate_heftcamb_rph.py",
    "planck_wrapper_sha256": "scripts/evaluate_planck_2018_fixed.py",
    "late_time_evaluator_sha256": "scripts/evaluate_kgb_late_time.py",
}
SUMMARY_PROVENANCE_KEYS = frozenset(
    (
        "generator_sha256",
        "planck_wrapper_sha256",
        "late_time_evaluator_sha256",
        "heftcamb_binary_sha256",
        "heftcamb_template_sha256",
    )
)
CACHE_INPUT_KEYS = (
    "alpha",
    "omega_m0",
    "omega_r0",
    "ombh2",
    "scalar_amp",
    "scalar_spectral_index",
    "tau",
    "a_planck",
)


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"expected a JSON object: {path}")
    return value


def exact_contract(config: Path) -> tuple[Path, dict[str, Any]]:
    loaded = yaml.safe_load(config.read_text(encoding="utf-8"))
    try:
        options = loaded["likelihood"][LIKELIHOOD_NAME]
    except (KeyError, TypeError) as error:
        raise RuntimeError(f"missing KGB likelihood definition: {config}") from error
    if not isinstance(options, dict):
        raise RuntimeError("KGB likelihood options are not a mapping")
    source_root = Path(options["source_root"]).resolve()
    components = {
        key: source_root / relative for key, relative in PROVENANCE_COMPONENTS.items()
    }
    components.update(
        {
            "joint_evaluator_sha256": source_root / "scripts/evaluate_kgb_joint_point.py",
            "heftcamb_binary_sha256": Path(options["heftcamb_binary"]).resolve(),
            "heftcamb_template_sha256": Path(options["heftcamb_template"]).resolve(),
        }
    )
    absent = [key for key, path in components.items() if not path.is_file()]
    if absent:
        raise FileNotFoundError("missing contract component: " + ", ".join(absent))
    contract = {
        "format": 1,
        **{key: sha256(path) for key, path in components.items()},
        "clipy_source": str(Path(options["clipy_source"]).resolve()),
        "planck_base": str(Path(options["planck_base"]).resolve()),
        "data_root": str(Path(options["data_root"]).resolve()),
        "spline_points": int(options["spline_points"]),
        "late_time_integration_nodes": int(options["late_time_integration_nodes"]),
    }
    return Path(options["output_root"]).resolve(), contract


def point_key(inputs: dict[str, Any], contract: dict[str, Any]) -> str:
    payload = json.dumps(
        {"inputs": inputs, "execution_contract": contract},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()[:24]


def audit_point(
    directory: Path, contract: dict[str, Any], *, deep_archive_verification: bool
) -> dict[str, str]:
    summary_path = directory / "joint_point.json"
    summary = read_json(summary_path)
    if summary.get("status") != "completed":
        raise RuntimeError(f"point is not completed: {directory}")
    inputs = summary.get("inputs")
    if not isinstance(inputs, dict):
        raise RuntimeError(f"point lacks immutable inputs: {directory}")
    cache_inputs = {name: inputs[name] for name in CACHE_INPUT_KEYS}
    if point_key(cache_inputs, contract) != directory.name:
        raise RuntimeError(f"cache key does not match inputs and execution contract: {directory}")
    provenance = summary.get("execution_provenance")
    if not isinstance(provenance, dict) or set(provenance) != SUMMARY_PROVENANCE_KEYS:
        raise RuntimeError(f"point lacks complete execution provenance: {directory}")
    expected_provenance = {
        key: contract[key] for key in SUMMARY_PROVENANCE_KEYS
    }
    if provenance != expected_provenance:
        raise RuntimeError(f"point provenance differs from declared execution contract: {directory}")
    likelihood = summary.get("minus_2_log_likelihood")
    if not isinstance(likelihood, dict):
        raise RuntimeError(f"point lacks likelihood components: {directory}")
    planck = float(likelihood["planck"])
    late_time = float(likelihood["late_time"])
    total = float(likelihood["total"])
    if not all(math.isfinite(value) for value in (planck, late_time, total)):
        raise RuntimeError(f"point has a non-finite likelihood: {directory}")
    if not math.isclose(total, planck + late_time, rel_tol=0.0, abs_tol=1.0e-10):
        raise RuntimeError(f"point likelihood components do not sum: {directory}")
    manifest_path = directory / MANIFEST
    archive_path = directory / ARCHIVE
    if not manifest_path.is_file() or not archive_path.is_file():
        raise RuntimeError(f"point is missing its verified artifact archive: {directory}")
    manifest = read_json(manifest_path)
    if sha256(archive_path) != manifest.get("archive_sha256"):
        raise RuntimeError(f"archive digest does not match manifest: {directory}")
    entries = manifest.get("files")
    if not isinstance(entries, list) or not entries:
        raise RuntimeError(f"archive manifest lacks file entries: {directory}")
    # The compactor verifies every member against this manifest before it
    # atomically writes the archive SHA-256. Rechecking that digest is enough
    # for a full production cache audit; ``--deep-archive-verification`` is
    # available for an independent member-by-member re-read on smaller sets.
    if deep_archive_verification:
        verify_archive(archive_path, entries)
    return {
        "point": directory.name,
        "joint_point_sha256": sha256(summary_path),
        "manifest_sha256": sha256(manifest_path),
        "archive_sha256": sha256(archive_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--cache-directory", type=Path)
    parser.add_argument(
        "--deep-archive-verification",
        action="store_true",
        help="Re-read and hash every retained member in addition to the archive digest.",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    configured_cache, contract = exact_contract(args.config.resolve())
    cache = args.cache_directory.resolve() if args.cache_directory else configured_cache
    if cache != configured_cache:
        raise RuntimeError("cache-directory differs from the immutable likelihood contract")
    if not cache.is_dir():
        raise FileNotFoundError(f"cache directory unavailable: {cache}")
    records = [
        audit_point(
            directory,
            contract,
            deep_archive_verification=args.deep_archive_verification,
        )
        for directory in sorted(cache.iterdir())
        if directory.is_dir() and (directory / "joint_point.json").is_file()
    ]
    record_bytes = json.dumps(records, sort_keys=True, separators=(",", ":")).encode("ascii")
    output = {
        "status": "passed",
        "model": "R-Universe covariant KGB",
        "scope": "Completed exact-point cache integrity audit; in-progress directories are excluded.",
        "config": str(args.config.resolve()),
        "config_sha256": sha256(args.config.resolve()),
        "cache_directory": str(cache),
        "execution_contract": contract,
        "archive_verification": (
            "archive SHA-256 plus compaction-time member validation"
            if not args.deep_archive_verification
            else "archive SHA-256 plus independent member-by-member validation"
        ),
        "completed_points": len(records),
        "point_record_sha256": hashlib.sha256(record_bytes).hexdigest(),
        "points": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
