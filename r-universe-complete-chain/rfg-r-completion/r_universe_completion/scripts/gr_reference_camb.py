#!/usr/bin/env python3
"""Execute a pinned GR CAMB reference calculation for the CMB data interface.

This script is deliberately a *GR reference only*. CAMB does not contain the
RFG-R extended operator ``bar_m5 deltaR3 deltaK`` and therefore no result from
this file may be labelled an RFG-R prediction or used in an RFG-R likelihood.
It records a reproducible end-to-end check of the ordinary photon, baryon,
CDM, neutrino, recombination, lensing and matter-transfer infrastructure.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "generated" / "tables"


@dataclass(frozen=True)
class Planck2018GRReference:
    """Published Planck 2018 base-LambdaCDM reference inputs.

    The first six entries are the TT,TE,EE+lowE+lensing posterior means in
    Planck Collaboration VI, A&A 641, A6 (2020), Table 2. The 0.06 eV
    neutrino mass is the fixed standard-mass assumption for the reference
    CAMB run. They define a regression calculation, not an RFG-R fit.
    """

    H0_km_s_Mpc: float = 67.36
    omega_b_h2: float = 0.02237
    omega_cdm_h2: float = 0.1200
    tau_reio: float = 0.0544
    ln_1e10_As: float = 3.044
    n_s: float = 0.9649
    sum_mnu_eV: float = 0.06
    T_CMB_K: float = 2.7255
    camb_nnu: float = 3.044

    @property
    def A_s(self) -> float:
        return math.exp(self.ln_1e10_As) * 1.0e-10


def _camb_module() -> Any:
    try:
        import camb
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "CAMB is required for the GR reference only. Install the pinned runtime with "
            "`python -m pip install camb==2.0.1`, then rerun this script."
        ) from error
    if camb.__version__ != "2.0.1":
        raise RuntimeError(f"expected CAMB 2.0.1, found {camb.__version__}")
    return camb


def calculate_reference(reference: Planck2018GRReference = Planck2018GRReference()) -> dict[str, object]:
    """Return directly computed GR spectra and transfer samples in CAMB units."""
    camb = _camb_module()
    from camb.bbn import get_predictor

    # The Planck likelihood preset uses this same public BBN interpolation.
    # The second argument is Delta N_eff, rather than N_eff itself.
    YHe_BBN = float(get_predictor().Y_He(reference.omega_b_h2, 0.0))
    parameters = camb.CAMBparams()
    parameters.set_cosmology(
        H0=reference.H0_km_s_Mpc,
        ombh2=reference.omega_b_h2,
        omch2=reference.omega_cdm_h2,
        mnu=reference.sum_mnu_eV,
        omk=0.0,
        tau=reference.tau_reio,
        TCMB=reference.T_CMB_K,
        nnu=reference.camb_nnu,
        YHe=YHe_BBN,
    )
    parameters.InitPower.set_params(As=reference.A_s, ns=reference.n_s)
    parameters.set_matter_power(redshifts=[0.0], kmax=2.0)
    parameters.set_for_lmax(2500, lens_potential_accuracy=1)
    results = camb.get_results(parameters)
    total = results.get_cmb_power_spectra(parameters, CMB_unit="muK")["total"]
    lensing = results.get_lens_potential_cls(lmax=2500)
    kh, redshifts, matter_power = results.get_matter_power_spectrum(minkh=0.01, maxkh=0.2, npoints=3)
    multipoles = (2, 30, 200, 800, 1200, 2000, 2500)
    return {
        "reference": reference,
        "camb_version": camb.__version__,
        "cmb": {ell: tuple(float(value) for value in total[ell, (0, 1, 3)]) for ell in multipoles},
        "lensing_pp": {ell: float(lensing[ell, 0]) for ell in multipoles},
        "matter": {float(k): float(matter_power[0, index]) for index, k in enumerate(kh)},
        "sigma8": float(results.get_sigma8_0()),
        "H0_check": float(results.hubble_parameter(0.0)),
        "YHe_BBN": YHe_BBN,
        "redshifts": tuple(float(redshift) for redshift in redshifts),
    }


def write_reference(result: dict[str, object], directory: Path = TABLES) -> tuple[Path, Path]:
    """Write machine-readable samples and a human-readable provenance record."""
    directory.mkdir(parents=True, exist_ok=True)
    reference = result["reference"]
    if not isinstance(reference, Planck2018GRReference):
        raise TypeError("reference metadata has the wrong type")
    rows: list[dict[str, object]] = []
    for ell, values in result["cmb"].items():
        tt, ee, te = values
        rows.append(
            {
                "observable": "D_ell_CMB",
                "coordinate": ell,
                "TT_uK2": tt,
                "EE_uK2": ee,
                "TE_uK2": te,
                "phi_phi_dimensionless": result["lensing_pp"][ell],
                "P_k_h3_Mpc3": "",
            }
        )
    for kh, value in result["matter"].items():
        rows.append(
            {
                "observable": "P_k_z0",
                "coordinate": kh,
                "TT_uK2": "",
                "EE_uK2": "",
                "TE_uK2": "",
                "phi_phi_dimensionless": "",
                "P_k_h3_Mpc3": value,
            }
        )
    csv_path = directory / "gr_camb_reference.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary_path = directory / "gr_camb_reference.md"
    metadata = asdict(reference)
    summary = [
        "# GR CAMB Reference Calculation",
        "",
        "This file records a numerical GR regression calculation only.",
        "It is not an RFG-R spectrum, likelihood, fit, or model comparison.",
        "",
        f"CAMB version: {result['camb_version']}",
        f"H0 check [km s^-1 Mpc^-1]: {result['H0_check']:.12e}",
        f"sigma8(z=0): {result['sigma8']:.12e}",
        f"BBN-consistent YHe: {result['YHe_BBN']:.15e}",
        "",
        "## Inputs",
        "",
    ]
    for key, value in metadata.items():
        summary.append(f"- {key}: {value}")
    summary.extend(
        [
            "",
            "The CMB rows are CAMB total lensed D_ell values in microkelvin squared.",
            "The lensing column uses CAMB's dimensionless lens-potential convention.",
            "The matter rows are linear P(k,z=0) in (Mpc/h)^3.",
            "YHe is evaluated by CAMB's public BBN_table_interpolator at the supplied",
            "omega_b h^2 and Delta N_eff=0, matching the Cobaya Planck reference preset.",
            "",
            "RFG-R requires a separate implementation retaining bar_m5 deltaR3 deltaK.",
        ]
    )
    summary_path.write_text("\n".join(summary) + "\n", encoding="utf-8")
    return csv_path, summary_path


def main() -> None:
    result = calculate_reference()
    csv_path, summary_path = write_reference(result)
    if not math.isclose(result["H0_check"], Planck2018GRReference().H0_km_s_Mpc, rel_tol=0.0, abs_tol=1.0e-12):
        raise RuntimeError("CAMB did not preserve the supplied H0")
    if not math.isfinite(result["sigma8"]) or result["sigma8"] <= 0.0:
        raise RuntimeError("nonphysical GR reference sigma8")
    print(f"Generated {csv_path}")
    print(f"Generated {summary_path}")
    print(f"CAMB {result['camb_version']} GR sigma8 = {result['sigma8']:.12e}")


if __name__ == "__main__":
    main()
