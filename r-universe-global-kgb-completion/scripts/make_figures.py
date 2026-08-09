#!/usr/bin/env python3
"""Generate a compact stability figure for the global R-alpha KGB completion."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ru_kgb import RUKGBParams, trajectory


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    rows = trajectory(RUKGBParams(), points=1801)
    a = np.array([row["a"] for row in rows])
    e = np.array([row["E"] for row in rows])
    omega_r = np.array([row["Omega_R"] for row in rows])
    b = np.array([row["alpha_B"] for row in rows])
    k = np.array([row["alpha_K"] for row in rows])
    qs = np.array([row["Q_s_over_Mpl2"] for row in rows])
    cs = np.array([row["c_s2"] for row in rows])

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.45))
    axes[0].loglog(a, e, color="#0f766e", lw=2.2, label=r"$E(a)$")
    axes[0].loglog(a, omega_r, color="#c2410c", lw=2.2, label=r"$\Omega_R(a)$")
    axes[0].axvline(1.0, color="0.4", lw=1.0, ls="--")
    axes[0].set_xlabel(r"scale factor $a$")
    axes[0].set_ylabel("background quantity")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].semilogx(a, b, color="#2563eb", lw=2.1, label=r"$\alpha_B$")
    axes[1].semilogx(a, k, color="#be123c", lw=2.1, label=r"$\alpha_K$")
    axes[1].semilogx(a, qs, color="#7c3aed", lw=2.1, label=r"$Q_s/M_{\rm Pl}^2$")
    axes[1].axvline(1.0, color="0.4", lw=1.0, ls="--")
    axes[1].set_xlabel(r"scale factor $a$")
    axes[1].set_ylabel("scalar stability coefficient")
    axes[1].legend(frameon=False, fontsize=8)

    axes[2].semilogx(a, cs, color="#15803d", lw=2.3, label=r"$c_s^2$")
    axes[2].axhline(1.0, color="0.4", lw=1.0, ls="--")
    axes[2].axvline(1.0, color="0.4", lw=1.0, ls="--")
    axes[2].set_ylim(0.96, 1.04)
    axes[2].set_xlabel(r"scale factor $a$")
    axes[2].set_ylabel(r"scalar sound speed squared")
    axes[2].legend(frameon=False, fontsize=8)

    for ax in axes:
        ax.grid(alpha=0.22)
    fig.tight_layout()
    output = ROOT / "generated" / "ru_kgb_stability.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={"Creator": "R-Universe KGB reproducibility package", "CreationDate": None, "ModDate": None})
    fig.savefig(output.with_suffix(".png"), dpi=200, metadata={"Software": "R-Universe KGB reproducibility package"})
    plt.close(fig)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
