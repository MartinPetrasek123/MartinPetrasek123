# Directly Testable RFG-RXi Factors

This calculation evaluates only quantities that follow directly from the
defined RFG-RXi action and the already stated local-GR matching rule. It is
not a CMB, matter-power, BAO, supernova, or joint cosmological likelihood.

For the declared unfitted benchmarks `Xi=1` and `Xi=2`, the script evaluates
the exact R-Universe background, the completed tensor coefficient

```text
Q_T^Xi(a) = Q_epsilon(a) + Xi,
c_T = 1,
dL_GW/dL_EM = sqrt(Q_T^Xi(1)/Q_T^Xi(a)),
```

and the local PPN sector. The Weyl switch is exactly zero in the Solar-System
matching domain, so both benchmarks have `gamma=beta=1` and
`alpha1=alpha2=0`. The Gaussian Cassini factor is therefore a valid data
factor for this action family:

```text
-2 ln L_Cassini = 0.8336483932.
```

The result is independent of Xi because Xi is multiplied by the same local
constant-zero Weyl switch. The DAE columns in the summary are read from the
separate exact action audit, which finds no root on the declared 49-by-81
grid for either benchmark.

Run, after the DAE audit:

```bash
cd r_universe_completion
python3 scripts/rfg_xi_completion.py
python3 scripts/rfg_xi_observables.py
python3 scripts/validate_rfg_xi_observables.py
```

The output files are `generated/tables/rfg_xi_observables.csv` and
`generated/tables/rfg_xi_validated_factors.csv`. A CMB likelihood would need
the still-unimplemented action-faithful photon--baryon--CDM--neutrino solver;
the script intentionally does not substitute a GR or LCDM transfer function.
