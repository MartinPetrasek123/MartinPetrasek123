#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "code" / "boltzmann_camb_results.json"
OUT_CLS = ROOT / "tables" / "camb_lensed_cls_lcdm_reference.csv"
OUT_PK = ROOT / "tables" / "camb_matter_power_lcdm_reference.csv"


def main():
    import camb
    from camb import model

    fits = json.loads((ROOT / "code" / "extended_results.json").read_text())
    lcdm = fits[1]["models"]["lcdm"]

    # Physical densities are Planck-like reference values. The late-time fit
    # supplies H0 and Omega_m; CAMB then computes the early-universe sound
    # horizon self-consistently for this reference sector.
    h = lcdm["H0"] / 100.0
    ombh2 = 0.02237
    omch2 = max(lcdm["Omega_m0"] * h * h - ombh2 - 0.00064, 1e-4)

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=lcdm["H0"], ombh2=ombh2, omch2=omch2, mnu=0.06, omk=0.0, tau=0.0544)
    pars.InitPower.set_params(As=2.1e-9, ns=0.9649, r=0.0)
    pars.set_for_lmax(2500, lens_potential_accuracy=1)
    pars.set_matter_power(redshifts=[0.0, 0.5, 1.0, 2.0], kmax=2.0)

    results = camb.get_results(pars)
    derived = results.get_derived_params()
    cls = results.get_cmb_power_spectra(pars, CMB_unit="muK")
    lensed = cls["total"]
    ell = np.arange(lensed.shape[0])
    np.savetxt(
        OUT_CLS,
        np.column_stack([ell, lensed[:, 0], lensed[:, 1], lensed[:, 2], lensed[:, 3]]),
        delimiter=",",
        header="ell,TT,EE,BB,TE",
        comments="",
    )

    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=2.0, npoints=180)
    rows = []
    for zi, pki in zip(z, pk):
        for kval, pkval in zip(kh, pki):
            rows.append((zi, kval, pkval))
    np.savetxt(
        OUT_PK,
        np.array(rows),
        delimiter=",",
        header="z,k_h_per_Mpc,P_k_Mpc_over_h_cubed",
        comments="",
    )

    out = {
        "engine": "CAMB",
        "camb_version": camb.__version__,
        "model": "LCDM reference early-universe calculation",
        "note": "This is not a full R1 Boltzmann implementation; it computes the reference early-universe sector and the fixed-rd replacement target.",
        "input": {
            "H0": lcdm["H0"],
            "Omega_m0": lcdm["Omega_m0"],
            "ombh2": ombh2,
            "omch2": omch2,
            "mnu_eV": 0.06,
            "tau": 0.0544,
            "As": 2.1e-9,
            "ns": 0.9649,
        },
        "derived": {k: float(v) for k, v in derived.items()},
        "outputs": {
            "lensed_cls": str(OUT_CLS.relative_to(ROOT)),
            "matter_power": str(OUT_PK.relative_to(ROOT)),
        },
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
