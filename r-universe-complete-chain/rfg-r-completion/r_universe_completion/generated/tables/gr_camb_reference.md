# GR CAMB Reference Calculation

This file records a numerical GR regression calculation only.
It is not an RFG-R spectrum, likelihood, fit, or model comparison.

CAMB version: 2.0.1
H0 check [km s^-1 Mpc^-1]: 6.736000000000e+01
sigma8(z=0): 8.110325278646e-01
BBN-consistent YHe: 2.456560393606866e-01

## Inputs

- H0_km_s_Mpc: 67.36
- omega_b_h2: 0.02237
- omega_cdm_h2: 0.12
- tau_reio: 0.0544
- ln_1e10_As: 3.044
- n_s: 0.9649
- sum_mnu_eV: 0.06
- T_CMB_K: 2.7255
- camb_nnu: 3.044

The CMB rows are CAMB total lensed D_ell values in microkelvin squared.
The lensing column uses CAMB's dimensionless lens-potential convention.
The matter rows are linear P(k,z=0) in (Mpc/h)^3.
YHe is evaluated by CAMB's public BBN_table_interpolator at the supplied
omega_b h^2 and Delta N_eff=0, matching the Cobaya Planck reference preset.

RFG-R requires a separate implementation retaining bar_m5 deltaR3 deltaK.
