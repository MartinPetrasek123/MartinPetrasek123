#!/usr/bin/env python3
"""Run one reproducible, independently seeded KGB Cobaya chain."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from cobaya.run import run
from cobaya.yaml import yaml_load_file


ROOT = Path(__file__).resolve().parents[1]
# Cobaya resolves the package-local likelihood by import name.  Make every
# direct or supervisor-launched invocation independent of shell PYTHONPATH.
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def configuration_path(config: Path, output: Path, resume: bool) -> Path:
    """Use the immutable input contract recorded for a resumed chain.

    Cobaya keeps the proposal adaptation state and checks the posterior model on
    resume.  Reading the original ``.input.yaml`` as well makes the contract
    explicit: later edits to the template cannot silently alter a historical
    chain.  ``max_samples`` remains an intentionally mutable stopping control.
    """
    saved_input = output.with_suffix(".input.yaml")
    if resume and saved_input.is_file():
        return saved_input
    return config


def configure_sampler_seed(info: dict, *, seed: int, resume: bool) -> None:
    """Bind a fresh chain to Cobaya's native RNG without changing a resume.

    Cobaya checkpoints proposal adaptation but not the bit-generator state.  A
    resumed chain therefore retains the exact sampler contract saved in its
    original input file; a new chain receives its declared independent seed.
    """
    sampler = info["sampler"]["mcmc"]
    recorded_seed = sampler.get("seed")
    if resume:
        # Production files created before this native seed binding have no
        # ``seed`` entry.  Preserve their saved contract exactly rather than
        # retroactively modifying it during an extension.
        if recorded_seed is not None and int(recorded_seed) != seed:
            raise RuntimeError(
                "resume seed differs from the immutable chain contract: "
                f"{recorded_seed} != {seed}"
            )
        return
    sampler["seed"] = seed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "inference" / "kgb_joint_planck_late.yaml")
    parser.add_argument("--output", type=Path, required=True, help="Cobaya output prefix for this independent chain")
    parser.add_argument(
        "--max-samples",
        type=int,
        required=True,
        help="Maximum stored Cobaya Markov states for this chain.",
    )
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--resume", action="store_true", help="Resume a compatible interrupted Cobaya chain")
    parser.add_argument("--force", action="store_true", help="Allow overwriting an existing output prefix")
    args = parser.parse_args()
    if args.max_samples < 2:
        raise ValueError("max-samples must be at least two")
    if args.resume and args.force:
        raise ValueError("resume and force cannot be used together")
    config_path = configuration_path(args.config, args.output, args.resume)
    if not config_path.is_file():
        raise FileNotFoundError(f"configuration unavailable: {config_path}")
    info = yaml_load_file(str(config_path))
    info["sampler"]["mcmc"]["max_samples"] = args.max_samples
    # Cobaya owns the proposal RNG.  Supplying its native seed option makes
    # each independent chain reproducible instead of only seeding NumPy's
    # legacy global generator, which Cobaya does not use for MCMC proposals.
    configure_sampler_seed(info, seed=args.seed, resume=args.resume)
    run(
        info,
        output=str(args.output),
        resume=args.resume,
        force=args.force,
        no_mpi=True,
    )


if __name__ == "__main__":
    main()
