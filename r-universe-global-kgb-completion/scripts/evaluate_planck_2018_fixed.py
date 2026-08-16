#!/usr/bin/env python3
"""Evaluate the official Planck 2018 likelihoods at one R-Universe point.

This is deliberately a *fixed-point* calculation.  The cosmological and
primordial parameters are those recorded in the H-EFTCAMB input file; it does
not optimize them and it must not be interpreted as a posterior or a model
comparison.  Its purpose is to exercise the full photon--baryon--CDM--neutrino
Boltzmann hierarchy and evaluate the resulting spectra with the distributed,
official Planck 2018 likelihood objects.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SPECTRA = ROOT / "generated" / "heftcamb"
PLANCK_ABSOLUTE_CALIBRATION_MEAN = 1.0
PLANCK_ABSOLUTE_CALIBRATION_SIGMA = 0.0025


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spectra-dir", type=Path, default=DEFAULT_SPECTRA,
        help="Directory containing the H-EFTCAMB output files.",
    )
    parser.add_argument(
        "--planck-base", type=Path,
        default=os.environ.get("PLANCK_2018_BASE"),
        help="Planck plc_3.0 directory containing the official likelihood files.",
    )
    parser.add_argument(
        "--clipy-source", type=Path,
        default=os.environ.get("CLIPY_SOURCE"),
        help="Directory containing the clipy package source.",
    )
    parser.add_argument(
        "--a-planck", type=float, default=1.0,
        help="Fixed Planck absolute-calibration nuisance parameter.",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "generated" / "planck_2018_fixed_loglike.json",
        help="JSON report path.",
    )
    args = parser.parse_args()
    if args.planck_base is None:
        parser.error("set PLANCK_2018_BASE or pass --planck-base")
    if args.clipy_source is None:
        parser.error("set CLIPY_SOURCE or pass --clipy-source")
    return args


def require_path(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def load_table(path: Path, columns: int) -> np.ndarray:
    require_path(path, "Spectrum file")
    values = np.loadtxt(path, comments="#")
    if values.ndim != 2 or values.shape[1] != columns:
        raise ValueError(f"{path} must have {columns} columns; found shape {values.shape}")
    ell = values[:, 0].astype(int)
    if not np.array_equal(ell, np.arange(ell[0], ell[-1] + 1)):
        raise ValueError(f"Multipoles in {path} are not contiguous")
    if ell[0] != 2:
        raise ValueError(f"{path} must begin at ell=2; begins at ell={ell[0]}")
    return values


def dl_to_cl(ell: np.ndarray, dl: np.ndarray) -> np.ndarray:
    """Convert CAMB's D_ell = ell(ell+1) C_ell/(2 pi) to raw C_ell."""
    return dl * (2.0 * np.pi) / (ell * (ell + 1.0))


def lensing_dl_to_cl(ell: np.ndarray, dl_phi_phi: np.ndarray) -> np.ndarray:
    """Convert CAMB's [ell(ell+1)]^2 C_phi_phi/(2 pi) to raw C_phi_phi."""
    return dl_phi_phi * (2.0 * np.pi) / (ell * (ell + 1.0)) ** 2


def cmb_cls(lensed: np.ndarray, required_lmax: int) -> np.ndarray:
    ell = lensed[:, 0].astype(int)
    if ell[-1] < required_lmax:
        raise ValueError(f"Lensed CMB output ends at {ell[-1]}, need {required_lmax}")
    cls = np.zeros((6, required_lmax + 1), dtype=np.float64)
    use = ell <= required_lmax
    index = ell[use]
    # Native CAMB file ordering: ell, D_TT, D_EE, D_BB, D_TE.
    cls[0, index] = dl_to_cl(index, lensed[use, 1])
    cls[1, index] = dl_to_cl(index, lensed[use, 2])
    cls[2, index] = dl_to_cl(index, lensed[use, 3])
    cls[3, index] = dl_to_cl(index, lensed[use, 4])
    return cls


def lensing_cls(unlensed_with_phi: np.ndarray, required_lmax: int) -> np.ndarray:
    ell = unlensed_with_phi[:, 0].astype(int)
    if ell[-1] < required_lmax:
        raise ValueError(f"Lens-potential output ends at {ell[-1]}, need {required_lmax}")
    cls = np.zeros((7, required_lmax + 1), dtype=np.float64)
    use = ell <= required_lmax
    index = ell[use]
    # Native CAMB lenspotential output: D_phi-phi=[ell(ell+1)]^2 C_phi-phi
    # /(2 pi), then D_TT, D_EE, D_BB, D_TE, D_Tphi and D_Ephi. Planck's
    # clik lensing order is phi-phi, TT, EE, BB, TE, TB, EB, all raw C_ell.
    cls[0, index] = lensing_dl_to_cl(index, unlensed_with_phi[use, 5])
    cls[1, index] = dl_to_cl(index, unlensed_with_phi[use, 1])
    cls[2, index] = dl_to_cl(index, unlensed_with_phi[use, 2])
    cls[3, index] = dl_to_cl(index, unlensed_with_phi[use, 3])
    cls[4, index] = dl_to_cl(index, unlensed_with_phi[use, 4])
    return cls


def scalar_value(value: Any) -> float:
    return float(np.asarray(value))


def absolute_calibration_prior(a_planck: float) -> dict[str, Any]:
    pull = (a_planck - PLANCK_ABSOLUTE_CALIBRATION_MEAN) / PLANCK_ABSOLUTE_CALIBRATION_SIGMA
    log_likelihood = -0.5 * pull * pull
    return {
        "log_likelihood": log_likelihood,
        "minus_2_log_likelihood": -2.0 * log_likelihood,
        "mean": PLANCK_ABSOLUTE_CALIBRATION_MEAN,
        "sigma": PLANCK_ABSOLUTE_CALIBRATION_SIGMA,
        "label": "Planck absolute-calibration Gaussian prior",
    }


def report_path(path: Path, external_label: str) -> str:
    """Keep public reports reproducible without exposing machine-specific paths."""
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return external_label


def evaluate(label: str, likelihood: Any, cls: np.ndarray, nuisance: dict[str, float]) -> dict[str, Any]:
    lmax = [int(item) for item in likelihood.get_lmax()]
    # Some official low-ell implementations multiply an entire input row by
    # their short ell grid, so their matrix must end exactly at its lmax.
    value = scalar_value(likelihood(cls[:, : max(lmax) + 1], nuisance))
    return {
        "log_likelihood": value,
        "minus_2_log_likelihood": -2.0 * value,
        "lmax": lmax,
        "extra_parameters": list(likelihood.get_extra_parameter_names()),
        "label": label,
    }


def main() -> None:
    args = parse_args()
    os.environ.setdefault("CLIPY_NOJAX", "1")
    require_path(args.clipy_source, "clipy source directory")
    require_path(args.planck_base, "Planck 2018 base directory")
    sys.path.insert(0, str(args.clipy_source))
    import clipy  # pylint: disable=import-outside-toplevel

    lensed = load_table(args.spectra_dir / "ru_kgb_rph_lensedCls.dat", 5)
    unlensed_with_phi = load_table(args.spectra_dir / "ru_kgb_rph_lenspotentialCls.dat", 8)
    nuisance = {"A_planck": args.a_planck}

    paths = {
        "plik_lite_TTTEEE": args.planck_base / "hi_l/plik_lite/plik_lite_v22_TTTEEE.clik",
        "commander_lowl_TT": args.planck_base / "low_l/commander/commander_dx12_v3_2_29.clik",
        "simall_lowl_EE": args.planck_base / "low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik",
        "lensing": args.planck_base / "lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8.clik_lensing",
    }
    for path in paths.values():
        require_path(path, "Planck likelihood")

    highl = clipy.clik(str(paths["plik_lite_TTTEEE"]))
    lowt = clipy.clik(str(paths["commander_lowl_TT"]))
    lowe = clipy.clik(str(paths["simall_lowl_EE"]))
    lens = clipy.clik_lensing(str(paths["lensing"]))

    cmb = cmb_cls(lensed, max(int(np.max(highl.get_lmax())), int(np.max(lowt.get_lmax())), int(np.max(lowe.get_lmax()))))
    lens_cls = lensing_cls(unlensed_with_phi, int(np.max(lens.get_lmax())))

    components = {
        "plik_lite_TTTEEE": evaluate("Planck 2018 Plik-lite TTTEEE", highl, cmb, nuisance),
        "commander_lowl_TT": evaluate("Planck 2018 Commander low-ell TT", lowt, cmb, nuisance),
        "simall_lowl_EE": evaluate("Planck 2018 SimAll low-ell EE", lowe, cmb, nuisance),
        "lensing": evaluate("Planck 2018 CMB-dependent lensing", lens, lens_cls, nuisance),
        "absolute_calibration_prior": absolute_calibration_prior(args.a_planck),
    }
    total_loglike = sum(item["log_likelihood"] for item in components.values())
    report = {
        "status": "completed",
        "scope": (
            "Official Planck 2018 likelihood evaluation at one fixed R-Universe "
            "KGB point, including the documented Planck absolute-calibration "
            "Gaussian prior A_planck=1.0000+-0.0025. This is not an optimization, "
            "posterior, evidence, or comparison to LambdaCDM."
        ),
        "inputs": {
            "spectra_dir": report_path(args.spectra_dir, "external spectrum directory supplied at runtime"),
            "planck_base": report_path(args.planck_base, "external Planck 2018 plc_3.0 distribution"),
            "a_planck_fixed": args.a_planck,
            "absolute_calibration_prior": {
                "mean": PLANCK_ABSOLUTE_CALIBRATION_MEAN,
                "sigma": PLANCK_ABSOLUTE_CALIBRATION_SIGMA,
            },
            "lensed_ell_max": int(lensed[-1, 0]),
            "lenspotential_ell_max": int(unlensed_with_phi[-1, 0]),
            "spectrum_convention": {
                "cmb": "CAMB D_ell converted to raw C_ell before clik evaluation",
                "phi_phi": "CAMB D_ell^phi_phi converted to raw C_ell^phi_phi before clik evaluation",
            },
        },
        "components": components,
        "total": {
            "log_likelihood": total_loglike,
            "minus_2_log_likelihood": -2.0 * total_loglike,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
