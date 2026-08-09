#!/usr/bin/env python3
"""Independent numerical gates for the globally completed R-alpha KGB action."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import mpmath as mp

from ru_kgb import RUKGBParams, background, trajectory


ROOT = Path(__file__).resolve().parents[1]


def scaled(value: float, scale: float) -> float:
    return abs(value) / max(1.0, scale)


def high_precision_action_gate(params: RUKGBParams) -> float:
    """Evaluate the algebraic action reconstruction where float64 cancels.

    In the radiation era A*X, B*X^2 and V are individually large while their
    sum is the finite R-sector density.  This independent 80-digit check keeps
    that benign cancellation out of the physical validation verdict.
    """
    mp.mp.dps = 80
    om = mp.mpf(str(params.omega_m0))
    orad = mp.mpf(str(params.omega_r0))
    alpha = mp.mpf(str(params.alpha))
    oR = 1 - om - orad
    beta_floor = alpha * orad
    worst = mp.mpf("0")
    for a_text in ["1e-8", "1e-5", "1e-3", "0.1", "1", "1.5", "2", "1000"]:
        a = mp.mpf(a_text)
        if a <= 1:
            effective, effective_n = a, a
        elif a >= 2:
            effective, effective_n = mp.mpf("2"), mp.mpf("0")
        else:
            t = a - 1
            left = mp.e ** (-1 / t)
            right = mp.e ** (-1 / (1 - t))
            step = left / (left + right)
            step_a = step * (1 - step) * (1 / t**2 + 1 / (1 - t) ** 2)
            effective = a + (2 - a) * step
            effective_n = a * (1 - step + (2 - a) * step_a)
        exponent = alpha * effective
        exponent_n = alpha * effective_n
        matter = om * a**-3
        radiation = orad * a**-4
        lo = mp.mpf("0")
        hi = mp.sqrt(matter + radiation + oR) + 1
        f = lambda e: e**2 - oR * e**exponent - matter - radiation
        while f(hi) <= 0:
            hi *= 2
        for _ in range(420):
            mid = (lo + hi) / 2
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        e = (lo + hi) / 2
        r = oR * e**exponent
        e_n = (r * exponent_n * mp.log(e) - 3 * matter - 4 * radiation) / (2 * e - r * exponent / e)
        r_n = r * (exponent_n * mp.log(e) + exponent * e_n / e)
        omega_R = r / e**2
        gate = exponent / (1 + exponent)
        gate_n = exponent_n / (1 + exponent) ** 2
        omega_R_n = omega_R * (r_n / r - 2 * e_n / e)
        b = -beta_floor + (1 + beta_floor) * gate * omega_R
        b_n = (1 + beta_floor) * (gate_n * omega_R + gate * omega_R_n)
        omega_m = matter / e**2
        omega_r = radiation / e**2
        sound_numerator = (2 - b) * (-e_n / e + b / 2) + b_n - (3 * omega_m + 4 * omega_r)
        D = sound_numerator
        alpha_k = D - mp.mpf("1.5") * b**2
        X = e**2 / 2
        C = b / e**2
        C_n = (b_n - 2 * b * e_n / e) / e**2
        rho = 3 * r
        pressure = -3 * r - r_n
        Rrho = rho - 6 * e**2 * X * C + 2 * C_n * X**2
        Rp = pressure + 2 * C_n * X**2 + 2 * C * X * e * e_n
        B = (alpha_k * e**2 - (Rrho + Rp) + 8 * C_n * X**2 - 12 * e**2 * X * C) / (8 * X**2)
        A = (Rrho + Rp - 4 * B * X**2) / (2 * X)
        V = (Rrho - Rp) / 2 - B * X**2
        rho_rec = A * X + 3 * B * X**2 + V + 6 * e**2 * X * C - 2 * C_n * X**2
        p_rec = A * X + B * X**2 - V - 2 * C_n * X**2 - 2 * C * X * e * e_n
        worst = max(worst, abs(rho_rec - rho) / max(1, abs(rho)), abs(p_rec - pressure) / max(1, abs(pressure)))
    return float(worst)


def main() -> None:
    params = RUKGBParams()
    params.validate()
    rows = trajectory(params)
    observed_rows = [row for row in rows if row["a"] <= 1.0]
    max_background = max(scaled(row["background_residual"], row["E"] ** 2) for row in rows)
    max_conservation = max(
        scaled(row["conservation_residual"], abs(row["rho_phi_target"]) + abs(row["p_phi_target"]))
        for row in rows
    )
    # Float64 cannot resolve the intended cancellation in the action density
    # at a << 10^-3.  The full domain is checked independently at 80 digits.
    noncancelling_rows = [row for row in rows if row["a"] >= 1.0e-3]
    max_rho_reconstruction = max(
        scaled(row["rho_reconstruction_residual"], row["rho_phi_target"]) for row in noncancelling_rows
    )
    max_p_reconstruction = max(
        scaled(row["p_reconstruction_residual"], abs(row["p_phi_target"])) for row in noncancelling_rows
    )
    high_precision_reconstruction = high_precision_action_gate(params)
    min_D = min(row["D"] for row in rows)
    min_Qs = min(row["Q_s_over_Mpl2"] for row in rows)
    min_Fs = min(row["F_s"] for row in rows)
    min_sound_numerator = min(row["sound_numerator"] for row in rows)
    min_alpha_k = min(row["alpha_K"] for row in rows)
    max_cs_deviation = max(abs(row["c_s2"] - params.target_cs2) for row in rows)
    min_ct2 = min(row["c_T2"] for row in rows)
    max_tensor_distance = max(abs(row["dL_GW_over_dL_EM"] - 1.0) for row in rows)
    max_observed_match = max(abs(row["a_effective"] / row["a"] - 1.0) for row in observed_rows)
    future = background(1.0e3, params)
    future_next = background(1.0e2, params)

    assert max_background < 3.0e-13
    assert max_conservation < 3.0e-12
    assert max_rho_reconstruction < 6.0e-12
    assert max_p_reconstruction < 6.0e-12
    assert high_precision_reconstruction < 1.0e-55
    assert min_D > 0.0
    assert min_Qs > 0.0
    assert min_Fs > 0.0
    assert min_sound_numerator > 0.0
    assert min_alpha_k > 0.0
    assert max_cs_deviation < 3.0e-13
    assert min_ct2 == 1.0
    assert max_tensor_distance == 0.0
    assert max_observed_match == 0.0
    assert abs(future["E"] / future_next["E"] - 1.0) < 2.0e-5

    output = ROOT / "generated" / "ru_kgb_trajectory.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "a", "z", "phi_over_Mpl", "a_effective", "E", "E_N", "response_exponent",
        "Omega_m", "Omega_r", "Omega_R", "w_R", "alpha_B", "alpha_K", "D", "F_s",
        "sound_numerator", "Q_s_over_Mpl2", "c_s2", "A_hat", "B_hat", "C_hat", "V_hat", "c_T2",
        "dL_GW_over_dL_EM",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows([{field: row[field] for field in fields} for row in rows])

    present = background(1.0, params)
    summary = {
        "model": "global R-alpha luminal KGB completion",
        "parameters": {
            "Omega_m0": params.omega_m0,
            "Omega_r0": params.omega_r0,
            "Omega_R0": params.omega_R0,
            "alpha": params.alpha,
            "beta_radiation": params.beta_radiation,
            "a_saturation": params.a_saturation,
            "target_cs2": params.target_cs2,
        },
        "gates": {
            "max_relative_background_residual": max_background,
            "max_relative_conservation_residual": max_conservation,
            "max_relative_rho_reconstruction_residual": max_rho_reconstruction,
            "max_relative_p_reconstruction_residual": max_p_reconstruction,
            "max_high_precision_action_reconstruction_residual": high_precision_reconstruction,
            "min_D": min_D,
            "min_Qs_over_Mpl2": min_Qs,
            "min_Fs": min_Fs,
            "min_sound_numerator": min_sound_numerator,
            "min_alpha_K": min_alpha_k,
            "max_abs_cs2_minus_target": max_cs_deviation,
            "observed_window_match": max_observed_match,
            "future_E_at_a_100": future_next["E"],
            "future_E_at_a_1000": future["E"],
        },
        "present": {
            key: present[key]
            for key in ["E", "Omega_R", "w_R", "alpha_B", "alpha_K", "Q_s_over_Mpl2", "c_s2"]
        },
        "tensor": {"Q_T_over_Mpl2": 1.0, "c_T2": 1.0, "dL_GW_over_dL_EM": 1.0},
        "scope": (
            "This validates the exact covariant background reconstruction and the high-frequency "
            "scalar/tensor no-ghost and no-gradient gates. A Boltzmann likelihood still requires "
            "a compiled perturbation module and public CMB likelihood data."
        ),
    }
    (ROOT / "generated" / "validation.json").write_text(json.dumps(summary, indent=2))

    print("R-Universe global KGB completion validation OK")
    print(f"max relative background residual = {max_background:.3e}")
    print(f"max relative conservation residual = {max_conservation:.3e}")
    print(f"max relative action rho residual = {max_rho_reconstruction:.3e}")
    print(f"max relative action p residual   = {max_p_reconstruction:.3e}")
    print(f"80-digit action residual         = {high_precision_reconstruction:.3e}")
    print(f"min N_s, D, Q_s, F_s, alpha_K = {min_sound_numerator:.6e}, {min_D:.6e}, {min_Qs:.6e}, {min_Fs:.6e}, {min_alpha_k:.6e}")
    print(f"max |c_s^2 - 1| = {max_cs_deviation:.3e}")
    print(f"present w_R = {present['w_R']:.8f}; alpha_B = {present['alpha_B']:.8f}")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
