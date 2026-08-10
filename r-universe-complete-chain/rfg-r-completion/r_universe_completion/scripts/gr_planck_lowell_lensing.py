#!/usr/bin/env python3
"""Evaluate official Planck low-ell and lensing likelihoods for the GR reference.

This is an executed data-interface regression test, not an RFG-R likelihood.
The exact RFG-R backend does not yet exist because stock CAMB/EFTCAMB omits
the nonzero ``bar_m5 deltaR3 deltaK`` operator. High-ell Plik is intentionally
excluded here: it requires its full foreground/calibration nuisance vector,
which must be sampled or fixed from a documented joint fit, never guessed.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from gr_reference_camb import Planck2018GRReference


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "generated" / "tables"
COMPONENTS = ("planck_2018_lowl.TT", "planck_2018_lowl.EE", "planck_2018_lensing.clik")


def calculate_loglikes(packages_path: Path) -> dict[str, float]:
    """Run the official data components with a BBN-consistent GR CAMB model."""
    try:
        from cobaya.cosmo_input import create_input
        from cobaya.model import get_model
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "Cobaya and CAMB are required. Install `cobaya==3.6.2 camb==2.0.1` and the Planck data packages."
        ) from error
    if not packages_path.is_dir():
        raise RuntimeError(f"Cobaya packages path does not exist: {packages_path}")
    reference = Planck2018GRReference()
    info: dict[str, Any] = create_input(
        theory="camb",
        like_cmb="planck_2018",
        primordial="SFSR",
        geometry="flat",
        hubble="H",
        matter="omegab_h2, omegac_h2",
        neutrinos="one_heavy_planck",
        dark_energy="lambda",
        bbn="consistency",
        reionization="std",
    )
    info.pop("sampler", None)
    info["packages_path"] = str(packages_path)
    info["likelihood"] = {component: None for component in COMPONENTS}
    for name, value in {
        "logA": reference.ln_1e10_As,
        "ns": reference.n_s,
        "H0": reference.H0_km_s_Mpc,
        "ombh2": reference.omega_b_h2,
        "omch2": reference.omega_cdm_h2,
        "tau": reference.tau_reio,
        "A_planck": 1.0,
    }.items():
        info["params"][name] = {"value": value, "drop": name == "logA"}
    model = get_model(info)
    loglikes, _ = model.loglikes({})
    if len(loglikes) != len(COMPONENTS):
        raise RuntimeError("unexpected number of Planck likelihood components")
    return {component: float(value) for component, value in zip(COMPONENTS, loglikes)}


def write_result(loglikes: dict[str, float], directory: Path = TABLES) -> tuple[Path, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    csv_path = directory / "gr_planck_2018_lowell_lensing.csv"
    rows = [
        {"component": component, "loglike": value, "minus_2loglike": -2.0 * value}
        for component, value in loglikes.items()
    ]
    rows.append(
        {
            "component": "combined_lowell_lensing",
            "loglike": sum(loglikes.values()),
            "minus_2loglike": -2.0 * sum(loglikes.values()),
        }
    )
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary_path = directory / "gr_planck_2018_lowell_lensing.md"
    summary_path.write_text(
        "# GR Planck 2018 Low-Ell And Lensing Reference\n\n"
        "This is an executed likelihood evaluation of the official Planck low-T, "
        "low-E and lensing components at a fixed GR/CAMB reference point. It is "
        "not an RFG-R likelihood, posterior, best fit, or model comparison.\n\n"
        "The likelihood uses a BBN-consistent helium fraction from Cobaya's "
        "standard one-heavy-neutrino Planck preset and the defined nominal "
        "calibration A_planck=1. The high-ell Plik likelihood is installed but not "
        "evaluated here, because its foreground and calibration parameters have not "
        "been fitted.\n\n"
        f"Combined log likelihood: {sum(loglikes.values()):.12e}\n"
        f"Combined -2 log likelihood: {-2.0 * sum(loglikes.values()):.12e}\n",
        encoding="utf-8",
    )
    return csv_path, summary_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packages-path", type=Path, required=True, help="Cobaya external packages directory")
    arguments = parser.parse_args()
    loglikes = calculate_loglikes(arguments.packages_path)
    csv_path, summary_path = write_result(loglikes)
    print(f"Generated {csv_path}")
    print(f"Generated {summary_path}")
    for component, value in loglikes.items():
        print(f"{component}: loglike = {value:.12e}")


if __name__ == "__main__":
    main()
