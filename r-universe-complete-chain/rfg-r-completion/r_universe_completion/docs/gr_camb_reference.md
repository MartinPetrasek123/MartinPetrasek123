# GR CAMB Reference And RFG-R Boundary

`scripts/gr_reference_camb.py` executes an independently reproducible,
numerical GR reference calculation through CAMB 2.0.1. It exercises the
standard photon, baryon, CDM, neutrino, recombination, lensing, and
matter-transfer machinery with published Planck 2018 base-LambdaCDM inputs.
The generated tables are a regression record for that external infrastructure.

It is not an RFG-R calculation. CAMB's GR action has no
`bar_m5 deltaR3 deltaK` contribution. Calling its output an RFG-R CMB spectrum
would silently set a nonzero coefficient of the defined RFG-R action to zero.

## Pinned Reference Inputs

The reference fixes `H0=67.36 km s^-1 Mpc^-1`,
`omega_b h^2=0.02237`, `omega_c h^2=0.1200`, `tau=0.0544`,
`ln(10^10 A_s)=3.044`, and `n_s=0.9649`. These are the Planck 2018
TT,TE,EE+lowE+lensing posterior means from Table 2 of Planck Collaboration VI,
A&A 641, A6 (2020). `sum_mnu=0.06 eV` and `T_CMB=2.7255 K` are fixed standard
assumptions. The
pinned CAMB/Cobaya Planck preset passes `nnu=3.044`; this is the code's
standard neutrino-count convention, while the physical reference is commonly
quoted as `N_eff=3.046`. They are not fitted RFG-R parameters.

For exactly the same GR reference point as the Planck likelihood interface,
the generator calculates `YHe` from CAMB's public `BBN_table_interpolator` at
the stated `omega_b h^2` and `Delta N_eff=0`; it records the resulting value in
the generated provenance table. It is a fixed BBN-consistency calculation, not
a parameter estimated from these data.

Run it with a Python environment containing the pinned public wheel:

```bash
python -m pip install 'camb==2.0.1'
python scripts/gr_reference_camb.py
```

The generated CMB rows use CAMB's total lensed `D_ell` convention. The matter
rows are linear `P(k,z=0)`. Neither table is a Planck likelihood evaluation:
the official likelihood consumes full theory spectra and nuisance parameters,
not a few exported reference samples.

## Exact RFG-R Requirement

The RFG-R background/ADM map is already computed by
`scripts/extended_eft_mapping.py`, including the nonzero
`bar_m5=-M_Pl^2 Q_X/(3H0)`. An exact RFG-R backend must add that operator to
the scalar constraint, scalar evolution, spatial-traceless/anisotropic-stress,
initial-condition, and stability systems before producing spectra. The kinetic
photon and neutrino hierarchies must then feed their stress-energy moments back
into those metric equations. No stock CAMB or H-EFTCAMB run satisfies that
condition.

The pinned regression is checked with:

```bash
python scripts/validate_gr_reference_camb.py
```
