# Backend Capability Audit

## Audited Source

The locally downloaded `EFTCAMB/EFTCAMB` source at commit
`16d9c4e9f85751e30efd0a53b177941713078904` was inspected on 2026-08-10.
It provides the ordinary EFTCAMB functions `Omega` and `Gamma1` through
`Gamma6`, a full CAMB matter/radiation hierarchy, recombination and lensing.
The source cache and perturbation equations contain no `bar_m5`,
`deltaR3 deltaK`, or equivalent seventh extended-EFT coefficient.

That absence is decisive for RFG-R. The exact RFG-R ADM map is

```math
bar_m5 = - M_Pl^2 Q_X/(3 H0),
```

and it is nonzero on the cosmological branch. Assigning it zero is a different
action. Therefore a successful stock H-EFTCAMB build can only be a GR/standard
EFT infrastructure reference, not an RFG-R spectrum calculation.

## Required Extension Before An RFG-R Spectrum

1. Add `bar_m5(a)` and its analytic action derivatives to the EFT cache and
   model interface.
2. Vary the `bar_m5 deltaR3 deltaK` operator in the scalar lapse, shift,
   trace and spatial-traceless field equations. The latter must include the
   photon and neutrino anisotropic stresses.
3. Couple the exact photon temperature/polarization hierarchy and the
   collisionless massless/massive-neutrino equations to the metric source
   moments. A finite perfect-fluid closure is not permitted.
4. Derive regular adiabatic initial conditions for this specific constrained
   system and verify its GR limit using identical primordial inputs.
5. Add the full constraint/ghost/gradient gate at every time and wavenumber,
   then compare spectra across integration tolerances and hierarchy
   resolutions before evaluating any likelihood.

The first step cannot be replaced by a mapping to an existing `Gamma1`--
`Gamma6` input: those functions parameterize a different operator basis. The
second through fourth steps require the complete action variation. Until they
are independently derived and regression-tested, this package must refuse to
label spectra or likelihoods as RFG-R.
