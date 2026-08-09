# R-Universe RFG-R Completion

This package defines RFG-R: a regular, multiscale effective completion of the
R-Universe preferred-foliation branch. It is designed to answer a precise
question: what theory is to be tested when one asks for matter growth, CMB
spectra, and Solar-System consistency?

RFG-R is not a rhetorical strengthening of the original branch. It specifies:

- one smooth nonlinear action on the cosmological branch;
- a reconstructed potential that preserves the R-Universe background to
  `O((epsilon/X)^p)`;
- a local GR matching domain that makes the PPN prediction unambiguous;
- universal minimal coupling of matter to the Jordan metric;
- a full linear Einstein-Boltzmann likelihood protocol, with stability checks
  before every likelihood evaluation.

The resulting physical picture is simple: a relational response is active on
homogeneous Hubble-scale geometry, while the locally resolved weak-field EFT is
General Relativity. The crossover is explicit rather than implicit.

## Model Definition

For `p=4`, `0<theta<2`, and `epsilon>0`, define

```math
R_epsilon(X)=Omega_R0 X^(p+2)/(X^p+epsilon^p)^(1+theta/p),
```

```math
Q_epsilon(X)=1-[Omega_R0/(1+theta)]X^p/(X^p+epsilon^p)^(1+theta/p).
```

The potential is fixed, not guessed:

```math
V_epsilon(X)=-3X integral_0^X ds/s^2 [s^2-R_epsilon(s)-s^2Q_epsilon(s)-s^3Q_epsilon'(s)].
```

The cosmological action is the RFG action with `Q,V` replaced by these
functions. The background equation is

```math
E^2=Omega_m0 a^-3+Omega_r0 a^-4+R_epsilon(E).
```

The local matching action is exactly Einstein-Hilbert plus minimally coupled
matter whenever the dimensionless Weyl indicator

```math
W=sqrt(abs(C_abcd C^abcd))/(H0/c)^2
```

exceeds `W_LOCAL_GR=1e9`. A C-infinity switching function is given in
`docs/completion_derivation.md`. FLRW has `W=0`; at one astronomical unit from
the Sun the code obtains `W about 6e22`.

## Contents

- `paper/R_Universe_RFG_R_Completion.pdf` - standalone paper.
- `docs/completion_derivation.md` - full mathematical definition and limits.
- `docs/likelihood_pipeline.md` - matter, CMB, and PPN likelihood protocol.
- `scripts/rfg_regularized.py` - background solver, potential reconstruction,
  and ADM coefficient table.
- `scripts/validate_completion.py` - independent numerical checks.
- `scripts/ppn_likelihood.py` - local GR matching and Cassini likelihood.
- `generated/` - tables and figures generated from the scripts.

## Reproduce

```bash
bash scripts/run_all.sh
```

The command checks the exact local limit, high-X recovery, potential
reconstruction, background closure, positive tensor normalization, and the PPN
domain. It then generates the tables and figures.

## Empirical Rule

The included CMB/matter protocol is a full-spectrum likelihood definition, not
a compressed-distance surrogate. It must be evaluated with the official Planck
likelihood data and a full 3+1 Einstein-Boltzmann backend before any statement
about empirical preference over LCDM is made. The package never substitutes a
background fit for a CMB fit.
