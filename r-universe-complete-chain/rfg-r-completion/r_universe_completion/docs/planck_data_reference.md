# Executed Planck Data Reference

The official Planck 2018 low-temperature, low-polarization, high-ell Plik and
lensing packages were installed through Cobaya 3.6.2 on 2026-08-10 under a
separate external-data directory. The package contains the public Planck data
and `clipy` likelihood implementation; it is not committed to this repository.

`scripts/gr_planck_lowell_lensing.py` has executed the low-T, low-E and
lensing likelihoods for the pinned GR/CAMB reference point. It uses the Planck
one-heavy-neutrino CAMB preset, a BBN-consistent helium prediction, and the
defined nominal calibration `A_planck=1`. Its CSV contains the individual and
combined log likelihood values.

The high-ell Plik package is deliberately not evaluated at arbitrary nuisance
values. It requires calibration, foreground, point-source, SZ and CIB
parameters. A valid high-ell value therefore belongs in a documented joint
fit. Picking defaults after seeing the cosmological input would be an
uncontrolled approximation, not a data calculation.

This GR result proves only that the external data interface runs. It does not
test RFG-R. RFG-R spectra and a corresponding Planck likelihood remain
forbidden until a solver retaining `bar_m5 deltaR3 deltaK` is implemented and
validated against the complete action variation.

Re-run the recorded data-interface regression with:

```bash
python scripts/validate_gr_planck_lowell_lensing.py --packages-path /path/to/cobaya-planck-2018
```
