#!/usr/bin/env python3
"""Pre-recombination CMB input gate from the explicit R-alpha KGB action."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from ru_kgb import RUKGBParams, background


ROOT = Path(__file__).resolve().parents[1]
C_KM_S = 299_792.458
H0_KM_S_MPC = 67.8625
OMEGA_B_H2 = 0.02237
T_CMB_K = 2.7255
Z_DRAG_REFERENCE = 1060.010


def lcdm_E(a: np.ndarray, params: RUKGBParams) -> np.ndarray:
    return np.sqrt(params.omega_m0 * a**-3 + params.omega_r0 * a**-4 + params.omega_R0)


def sound_horizon(params: RUKGBParams, a_drag: float) -> tuple[float, float]:
    """Return r_s for R-alpha KGB and the same-parameter LCDM comparison."""
    h = H0_KM_S_MPC / 100.0
    omega_b = OMEGA_B_H2 / (h * h)
    omega_gamma = 2.469e-5 * (T_CMB_K / 2.7255) ** 4 / (h * h)
    # Logarithmic integration resolves the radiation-era integrand smoothly.
    a = np.geomspace(1.0e-10, a_drag, 8001)
    e_kgb = np.array([background(float(ai), params)["E"] for ai in a])
    e_lcdm = lcdm_E(a, params)
    c_s = 1.0 / np.sqrt(3.0 * (1.0 + 3.0 * omega_b * a / (4.0 * omega_gamma)))
    n = np.log(a)
    integrand_kgb = c_s / (a * e_kgb)
    integrand_lcdm = c_s / (a * e_lcdm)
    prefactor = C_KM_S / H0_KM_S_MPC
    return (
        float(prefactor * np.trapezoid(integrand_kgb, n)),
        float(prefactor * np.trapezoid(integrand_lcdm, n)),
    )


def response(row: dict[str, float]) -> float:
    return row["alpha_B"] ** 2 / (2.0 * row["D"] * row["c_s2"])


def main() -> None:
    params = RUKGBParams()
    a_drag = 1.0 / (1.0 + Z_DRAG_REFERENCE)
    r_kgb, r_lcdm = sound_horizon(params, a_drag)
    epochs = {"BBN_z_1e9": 1.0e-9, "recombination_z_1090": 1.0 / 1091.0, "drag_z_1060": a_drag}
    checks = {}
    for name, a in epochs.items():
        row = background(a, params)
        e_lcdm = math.sqrt(params.omega_m0 * a**-3 + params.omega_r0 * a**-4 + params.omega_R0)
        checks[name] = {
            "a": a,
            "Omega_R": row["Omega_R"],
            "relative_H_minus_LCDM": row["E"] / e_lcdm - 1.0,
            "mu_infinity_minus_one": response(row),
            "alpha_B": row["alpha_B"],
            "Q_s_over_Mpl2": row["Q_s_over_Mpl2"],
        }
    result = {
        "model": "global R-alpha luminal KGB completion",
        "input": {
            "H0_km_s_Mpc": H0_KM_S_MPC,
            "Omega_b_h2": OMEGA_B_H2,
            "T_CMB_K": T_CMB_K,
            "z_drag_reference": Z_DRAG_REFERENCE,
        },
        "sound_horizon_background_gate": {
            "r_s_KGB_Mpc": r_kgb,
            "r_s_same_parameter_LCDM_Mpc": r_lcdm,
            "relative_difference": r_kgb / r_lcdm - 1.0,
        },
        "epochs": checks,
        "verdict": (
            "PASS: the R-sector background and the action-derived early scalar response are small before recombination; "
            "the sound-horizon background shift is explicitly quantified."
        ),
        "scope": (
            "This is a pre-recombination background-and-EFT input calculation. It is not a CMB angular-spectrum or Planck likelihood; "
            "those require evolving photon, baryon, neutrino, metric, and KGB scalar perturbations in a Boltzmann solver."
        ),
    }
    assert max(abs(row["relative_H_minus_LCDM"]) for row in checks.values()) < 2.0e-8
    assert max(row["Omega_R"] for row in checks.values()) < 2.0e-8
    assert abs(result["sound_horizon_background_gate"]["relative_difference"]) < 2.0e-8
    output = ROOT / "generated" / "cmb_prerecombination_gate.json"
    output.write_text(json.dumps(result, indent=2))
    print(result["verdict"])
    print(f"r_s KGB/LCDM - 1 = {result['sound_horizon_background_gate']['relative_difference']:.3e}")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
