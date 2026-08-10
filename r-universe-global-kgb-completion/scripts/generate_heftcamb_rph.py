#!/usr/bin/env python3
"""Generate an H-EFTCAMB RPH input from the exact R-Universe KGB action.

The RPH interface evolves a Horndeski model specified by w_DE(a), alpha_K(a),
alpha_B(a), alpha_M(a), and alpha_T(a).  This script supplies the first three
from ru_kgb.py, with alpha_M=alpha_T=0 exactly.  The spline is only a numerical
representation for the external solver: the defining model remains ru_kgb.py.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from ru_kgb import RUKGBParams, background


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "generated" / "heftcamb"


def radiation_h(params: RUKGBParams, temperature: float = 2.7255, neff: float = 3.046) -> float:
    """Fix H0 from the action's Omega_r0 and a massless standard neutrino bath."""
    omega_gamma = 2.4728e-5 * (temperature / 2.7255) ** 4
    omega_r = omega_gamma * (1.0 + 0.22710731766 * neff)
    return math.sqrt(omega_r / params.omega_r0)


def spline_lines(prefix: str, a: np.ndarray, values: np.ndarray, null_value: float) -> list[str]:
    """Emit an RPH spline with its exterior value fixed by the KGB limit."""
    lines = [f"{prefix}_Spline_Pixels = {len(a)}", f"{prefix}_null_value = {null_value:.17e}"]
    lines.extend(f"{prefix}x{i} = {x:.17e}" for i, x in enumerate(a, start=1))
    lines.extend(f"{prefix}v{i} = {y:.17e}" for i, y in enumerate(values, start=1))
    return lines


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=401)
    parser.add_argument("--a-min", type=float, default=1.0e-10)
    parser.add_argument("--a-max", type=float, default=1.01)
    parser.add_argument("--grid", choices=("linear", "power", "log"), default="power")
    parser.add_argument("--power", type=float, default=3.0)
    parser.add_argument("--turn-on", type=float, default=1.0e-4)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.points < 6:
        raise ValueError("H-EFTCAMB spline requires at least six points")
    if not 0.0 < args.a_min < 1.0 < args.a_max:
        raise ValueError("require 0 < a-min < 1 < a-max so the solver endpoints lie inside the spline")
    if args.grid == "power" and args.power <= 1.0:
        raise ValueError("power grid exponent must exceed one")
    if not args.a_min < args.turn_on < 1.0:
        raise ValueError("require a-min < turn-on < 1")

    params = RUKGBParams()
    # RPH differentiates a cubic spline with respect to a itself.  A logarithmic
    # spacing therefore creates artificial derivatives at early times; use a
    # uniform a grid and calculate every node from the exact KGB construction.
    # The interval extends beyond the solver endpoints because RPH returns its
    # null value exactly at a spline endpoint.
    if args.grid == "linear":
        a = np.linspace(args.a_min, args.a_max, args.points)
    elif args.grid == "power":
        u = np.linspace(args.a_min ** (1.0 / args.power), args.a_max ** (1.0 / args.power), args.points)
        a = u**args.power
    else:
        a = np.geomspace(args.a_min, args.a_max, args.points)
    rows = [background(float(ai), params) for ai in a]
    w = np.array([row["w_R"] for row in rows])
    alpha_k = np.array([row["alpha_K"] for row in rows])
    alpha_b_bs = np.array([row["alpha_B"] for row in rows])
    # H-EFTCAMB's RPH variable uses alpha_B^EFTCAMB=-alpha_B^BS/2;
    # ru_kgb.py and the action validation use the Bellini--Sawicki convention.
    alpha_b_rph = -0.5 * alpha_b_bs
    h = radiation_h(params)
    ombh2 = 0.02237
    ommh2 = params.omega_m0 * h * h
    omch2 = ommh2 - ombh2
    if omch2 <= 0.0:
        raise RuntimeError("declared baryon density exceeds the action-defined matter density")

    args.output.mkdir(parents=True, exist_ok=True)
    function_table = args.output / "ru_kgb_rph_functions.csv"
    np.savetxt(
        function_table,
        np.column_stack([a, w, alpha_k, alpha_b_bs, alpha_b_rph]),
        delimiter=",",
        header="a,w_R,alpha_K,alpha_B_BelliniSawicki,RPH_alpha_B_EFTCAMB",
        comments="",
    )

    lines = [
        "# Generated from scripts/ru_kgb.py; do not edit numerical functions by hand.",
        "# RPH is the luminal Horndeski interface: alpha_M=alpha_T=0.",
        "output_file_headers = T",
        "EFTCAMB_write_background = T",
        "output_root = ru_kgb_rph",
        "get_scalar_cls = T",
        "get_vector_cls = F",
        "get_tensor_cls = F",
        "get_transfer = T",
        "do_lensing = T",
        "do_nonlinear = 0",
        "l_max_scalar = 2700",
        "k_eta_max_scalar = 10800",
        "use_physical = T",
        f"ombh2 = {ombh2:.17e}",
        f"omch2 = {omch2:.17e}",
        "omnuh2 = 0",
        "omk = 0",
        f"hubble = {100.0 * h:.17e}",
        "temp_cmb = 2.7255",
        "helium_fraction = 0.24",
        "massless_neutrinos = 3.046",
        "nu_mass_eigenstates = 0",
        "massive_neutrinos = 0",
        "initial_power_num = 1",
        "pivot_scalar = 0.05",
        "scalar_amp(1) = 2.1e-9",
        "scalar_spectral_index(1) = 0.9649",
        "scalar_nrun(1) = 0",
        "scalar_nrunrun(1) = 0",
        "reionization = T",
        "re_use_optical_depth = T",
        "re_optical_depth = 0.0544",
        "RECFAST_fudge = 1.14",
        "RECFAST_fudge_He = 0.86",
        "RECFAST_Heswitch = 6",
        "RECFAST_Hswitch = T",
        "initial_condition = 1",
        "CMB_outputscale = 7.42835025e12",
        "transfer_high_precision = T",
        "transfer_kmax = 2",
        "transfer_k_per_logint = 0",
        "transfer_num_redshifts = 1",
        "transfer_redshift(1) = 0",
        "transfer_filename(1) = transfer_out.dat",
        "transfer_matterpower(1) = matterpower.dat",
        "transfer_power_var = 7",
        "scalar_output_file = scalCls.dat",
        "lensed_output_file = lensedCls.dat",
        "lens_potential_output_file = lenspotentialCls.dat",
        "highL_unlensed_cl_template = HighLExtrapTemplate_lenspotentialCls.dat",
        "feedback_level = 1",
        "derived_parameters = T",
        "accurate_polarization = T",
        "accurate_reionization = T",
        "do_late_rad_truncation = F",
        "massive_nu_approx = 0",
        "EFTflag = 2",
        "AltParEFTmodel = 1",
        "RPHwDE = 9",
        "RPHintegratefromtoday = T",
        "RPHusealphaM = F",
        "RPHmassPmodel = 0",
        "RPHmassPmodel_ODE = 0",
        "RPHkineticitymodel = 9",
        "RPHkineticitymodel_ODE = 0",
        "RPHbraidingmodel = 9",
        "RPHbraidingmodel_ODE = 0",
        "RPHtensormodel = 0",
        "RPHtensormodel_ODE = 0",
        "EFT_ghost_math_stability = T",
        # EFTCAMB marks this auxiliary exponential-mode heuristic deprecated.
        # Physical Horndeski ghost and gradient gates remain enabled below.
        "EFT_mass_math_stability = F",
        "EFT_ghost_stability = T",
        "EFT_gradient_stability = T",
        "EFT_mass_stability = F",
        "EFT_additional_priors = T",
        # The exact R-sector coefficients are < 6e-13 at a=1e-4.  Earlier
        # than this their sign is below double precision in external EFT code;
        # the solver therefore evolves standard adiabatic initial conditions
        # until the exact KGB functions are numerically resolvable.
        f"EFTCAMB_turn_on_time = {args.turn_on:.17e}",
        f"EFTCAMB_stability_time = {args.turn_on:.17e}",
        "EFTCAMB_stability_threshold = 0.0",
        "model_background_num_points = 6000",
        "model_background_a_ini = 1.e-8",
        "model_background_a_final = 1.0",
        *spline_lines("RPHw", a, w, -1.0),
        *spline_lines("RPHkineticity", a, alpha_k, alpha_k[0]),
        *spline_lines("RPHbraiding", a, alpha_b_rph, alpha_b_rph[0]),
    ]
    ini = args.output / "ru_kgb_rph.ini"
    ini.write_text("\n".join(lines) + "\n", encoding="ascii")
    metadata = args.output / "ru_kgb_rph_metadata.txt"
    metadata.write_text(
        "\n".join(
            [
                "Exact source: scripts/ru_kgb.py",
                f"spline_points = {args.points}",
                f"a_min = {args.a_min:.17e}",
                f"a_max = {args.a_max:.17e}",
                f"grid = {args.grid}",
                f"grid_power = {args.power:.17e}",
                f"perturbation_turn_on_a = {args.turn_on:.17e}",
                f"H0_from_Omega_r_km_s_Mpc = {100.0 * h:.17e}",
                f"ombh2 = {ombh2:.17e}",
                f"omch2 = {omch2:.17e}",
                "neutrino_sector = three massless standard neutrinos",
                "braiding_convention = RPH alpha_B = -BelliniSawicki alpha_B / 2",
                "primordial_sector = fixed spectrum gate, not an inferred likelihood point",
            ]
        )
        + "\n",
        encoding="ascii",
    )
    print(ini)
    print(function_table)
    print(metadata)


if __name__ == "__main__":
    main()
