#!/usr/bin/env python3
"""Scan physical KGB stability discriminants around the R-alpha fit."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from ru_kgb import RUKGBParams, trajectory


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    # This deliberately covers the profile interval and extends beyond it.
    alphas = np.linspace(0.02, 0.95, 16)
    omegas_m = np.linspace(0.24, 0.38, 15)
    rows: list[dict[str, float]] = []
    for alpha in alphas:
        for omega_m in omegas_m:
            params = RUKGBParams(alpha=float(alpha), omega_m0=float(omega_m))
            samples = trajectory(params, points=801)
            min_N = min(samples, key=lambda row: row["sound_numerator"])
            min_D = min(samples, key=lambda row: row["D"])
            min_Qs = min(samples, key=lambda row: row["Q_s_over_Mpl2"])
            min_cs2 = min(samples, key=lambda row: row["c_s2"])
            rows.append({
                "alpha": float(alpha),
                "Omega_m0": float(omega_m),
                "min_sound_numerator": min_N["sound_numerator"],
                "a_at_min_sound_numerator": min_N["a"],
                "min_D": min_D["D"],
                "min_Qs_over_Mpl2": min_Qs["Q_s_over_Mpl2"],
                "a_at_min_Qs": min_Qs["a"],
                "min_cs2": min_cs2["c_s2"],
            })

    assert all(row["min_sound_numerator"] > 0.0 for row in rows)
    assert all(row["min_D"] > 0.0 for row in rows)
    assert all(row["min_Qs_over_Mpl2"] > 0.0 for row in rows)
    assert all(row["min_cs2"] > 0.0 for row in rows)
    output = ROOT / "generated" / "stability_scan.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "range": {
            "alpha": [float(alphas.min()), float(alphas.max())],
            "Omega_m0": [float(omegas_m.min()), float(omegas_m.max())],
            "a": [1.0e-8, 1.0e3],
            "samples": len(rows),
        },
        "global_minima": {
            "sound_numerator": min(row["min_sound_numerator"] for row in rows),
            "D": min(row["min_D"] for row in rows),
            "Qs_over_Mpl2": min(row["min_Qs_over_Mpl2"] for row in rows),
            "cs2": min(row["min_cs2"] for row in rows),
        },
        "verdict": "PASS: physical scalar stability discriminants are positive throughout the scan",
    }
    (ROOT / "generated" / "stability_scan.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
