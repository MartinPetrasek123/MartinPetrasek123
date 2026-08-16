#!/usr/bin/env python3
"""Summarize independent exact KGB Cobaya chains without altering samples."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import norm, rankdata


PARAMETERS = ("alpha", "omega_m0", "omega_r0", "ombh2", "logA", "ns", "tau", "A_planck")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chains", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--burn-in-fraction", type=float, default=0.25)
    parser.add_argument("--rhat-threshold", type=float, default=1.01)
    parser.add_argument("--minimum-bulk-ess", type=float, default=400.0)
    parser.add_argument("--minimum-tail-ess", type=float, default=400.0)
    return parser.parse_args()


def read_chain(path: Path, burn_fraction: float) -> tuple[dict[str, np.ndarray], int, int]:
    if not 0.0 <= burn_fraction < 1.0:
        raise ValueError("burn-in-fraction must lie in [0, 1)")
    lines = path.read_text(encoding="utf-8").splitlines()
    header = next(line for line in lines if line.startswith("#")).lstrip("#").split()
    rows = [line.split() for line in lines if line.strip() and not line.startswith("#")]
    if not rows:
        raise RuntimeError(f"chain has no samples: {path}")
    values = np.asarray(rows, dtype=float)
    data = {name: values[:, index] for index, name in enumerate(header)}
    weights = data["weight"].astype(int)
    total = int(weights.sum())
    cutoff = int(np.floor(total * burn_fraction))
    kept: dict[str, list[float]] = {name: [] for name in header}
    consumed = 0
    for row_index, weight in enumerate(weights):
        start = max(0, cutoff - consumed)
        retained = int(weight) - start
        if retained > 0:
            for name in header:
                kept[name].extend([float(data[name][row_index])] * retained)
        consumed += int(weight)
    return {name: np.asarray(values, dtype=float) for name, values in kept.items()}, total, total - cutoff


def quantile(values: np.ndarray, probability: float) -> float:
    return float(np.quantile(values, probability, method="linear"))


def split_chains(chains: list[np.ndarray]) -> np.ndarray | None:
    if len(chains) < 2:
        return None
    split: list[np.ndarray] = []
    for chain in chains:
        half = len(chain) // 2
        if half < 2:
            return None
        split.extend((chain[:half], chain[-half:]))
    count = min(len(chain) for chain in split)
    return np.asarray([chain[:count] for chain in split])


def rhat_from_split_samples(samples: np.ndarray | None) -> float | None:
    if samples is None:
        return None
    count = samples.shape[1]
    means = samples.mean(axis=1)
    within = samples.var(axis=1, ddof=1).mean()
    between = count * means.var(ddof=1)
    if within <= 0.0:
        return None
    return float(np.sqrt(((count - 1.0) / count * within + between / count) / within))


def rank_normalize(values: np.ndarray) -> np.ndarray:
    flattened = values.reshape(-1)
    ranks = rankdata(flattened, method="average")
    # Blom's offset avoids infinities at the two endpoints.
    probabilities = (ranks - 3.0 / 8.0) / (len(ranks) + 1.0 / 4.0)
    return norm.ppf(probabilities).reshape(values.shape)


def rank_normalized_split_rhat(chains: list[np.ndarray]) -> float | None:
    samples = split_chains(chains)
    if samples is None:
        return None
    rank_rhat = rhat_from_split_samples(rank_normalize(samples))
    folded = np.abs(samples - np.median(samples))
    folded_rhat = rhat_from_split_samples(rank_normalize(folded))
    values = [value for value in (rank_rhat, folded_rhat) if value is not None]
    return max(values) if values else None


def autocovariance(samples: np.ndarray) -> np.ndarray:
    """Return biased within-chain autocovariances using an FFT convolution."""
    count = samples.shape[1]
    centered = samples - samples.mean(axis=1, keepdims=True)
    transform = np.fft.rfft(centered, n=2 * count, axis=1)
    return np.fft.irfft(transform * transform.conjugate(), n=2 * count, axis=1)[:, :count] / count


def ess_from_split_samples(samples: np.ndarray | None) -> float | None:
    """Geyer's initial-positive, monotone sequence ESS for split chains."""
    if samples is None:
        return None
    chains, draws = samples.shape
    if draws < 4:
        return None
    within = samples.var(axis=1, ddof=1).mean()
    if not math.isfinite(within) or within <= 0.0:
        return None
    between = draws * samples.mean(axis=1).var(ddof=1)
    variance_plus = (draws - 1.0) / draws * within + between / draws
    if not math.isfinite(variance_plus) or variance_plus <= 0.0:
        return None
    mean_autocovariance = autocovariance(samples).mean(axis=0)
    rho = 1.0 - (within - mean_autocovariance) / variance_plus
    rho[0] = 1.0

    pair_sums: list[float] = []
    for start in range(0, draws - 1, 2):
        pair = float(rho[start] + rho[start + 1])
        if pair < 0.0:
            break
        if pair_sums:
            pair = min(pair, pair_sums[-1])
        pair_sums.append(pair)
    if not pair_sums:
        return float(chains * draws)
    integrated_autocorrelation = max(-1.0 + 2.0 * sum(pair_sums), 1.0)
    return float(min(chains * draws, chains * draws / integrated_autocorrelation))


def rank_normalized_bulk_ess(chains: list[np.ndarray]) -> float | None:
    samples = split_chains(chains)
    return ess_from_split_samples(rank_normalize(samples) if samples is not None else None)


def tail_ess(chains: list[np.ndarray]) -> float | None:
    samples = split_chains(chains)
    if samples is None:
        return None
    lower, upper = np.quantile(samples, (0.05, 0.95))
    lower_ess = ess_from_split_samples((samples <= lower).astype(float))
    upper_ess = ess_from_split_samples((samples >= upper).astype(float))
    values = [value for value in (lower_ess, upper_ess) if value is not None]
    return min(values) if values else None


def main() -> None:
    args = parse_args()
    if (
        not 0.0 <= args.burn_in_fraction < 1.0
        or not 1.0 < args.rhat_threshold <= 1.1
        or args.minimum_bulk_ess <= 0.0
        or args.minimum_tail_ess <= 0.0
    ):
        raise ValueError("invalid posterior diagnostic thresholds")
    parsed = [read_chain(path, args.burn_in_fraction) for path in args.chains]
    combined: dict[str, np.ndarray] = {}
    for name in PARAMETERS + ("chi2__likelihoods.kgb_joint_likelihood.KGBJointLikelihood",):
        combined[name] = np.concatenate([data[name] for data, _, _ in parsed])
    summaries: dict[str, Any] = {}
    for name in PARAMETERS:
        values = combined[name]
        summaries[name] = {
            "mean": float(values.mean()),
            "standard_deviation": float(values.std(ddof=1)),
            "median": quantile(values, 0.5),
            "lower_68_percent": quantile(values, 0.16),
            "upper_68_percent": quantile(values, 0.84),
            "split_rhat": rhat_from_split_samples(
                split_chains([data[name] for data, _, _ in parsed])
            ),
            "rank_normalized_split_rhat": rank_normalized_split_rhat(
                [data[name] for data, _, _ in parsed]
            ),
            "bulk_ess": rank_normalized_bulk_ess([data[name] for data, _, _ in parsed]),
            "tail_ess": tail_ess([data[name] for data, _, _ in parsed]),
        }
    chi2 = combined["chi2__likelihoods.kgb_joint_likelihood.KGBJointLikelihood"]
    output = {
        "model": "R-Universe covariant KGB",
        "scope": "Exact-chain summary after an explicitly declared per-chain burn-in removal. This is a posterior result only when the reported independent-chain convergence diagnostics are satisfactory.",
        "diagnostic_method": {
            "split_rhat": "rank-normalized folded split-Rhat",
            "bulk_ess": "rank-normalized split-chain ESS with Geyer's initial-positive monotone sequence",
            "tail_ess": "minimum of 5-percent and 95-percent indicator ESS with Geyer's initial-positive monotone sequence",
        },
        "publication_gate": {
            "rank_normalized_split_rhat_less_than": args.rhat_threshold,
            "minimum_bulk_ess": args.minimum_bulk_ess,
            "minimum_tail_ess": args.minimum_tail_ess,
        },
        "chains": [str(path) for path in args.chains],
        "burn_in_fraction": args.burn_in_fraction,
        "per_chain_total_weight": [total for _, total, _ in parsed],
        "per_chain_retained_weight": [retained for _, _, retained in parsed],
        "combined_retained_weight": int(len(chi2)),
        "minimum_minus_2_log_likelihood_in_retained_samples": float(chi2.min()),
        "maximum_rank_normalized_split_rhat": max(
            float(summary["rank_normalized_split_rhat"])
            for summary in summaries.values()
            if summary["rank_normalized_split_rhat"] is not None
        ),
        "minimum_bulk_ess": min(
            float(summary["bulk_ess"])
            for summary in summaries.values()
            if summary["bulk_ess"] is not None
        ),
        "minimum_tail_ess": min(
            float(summary["tail_ess"])
            for summary in summaries.values()
            if summary["tail_ess"] is not None
        ),
        "passes_publication_gate": all(
            summary["rank_normalized_split_rhat"] is not None
            and summary["bulk_ess"] is not None
            and summary["tail_ess"] is not None
            and float(summary["rank_normalized_split_rhat"]) < args.rhat_threshold
            and float(summary["bulk_ess"]) >= args.minimum_bulk_ess
            and float(summary["tail_ess"]) >= args.minimum_tail_ess
            for summary in summaries.values()
        ),
        "parameters": summaries,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
