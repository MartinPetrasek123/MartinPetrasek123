# R-Universe Locked Finite-Window Theory

Standalone finite-window calibration theory and internal viability gate.

## Main result

The hard Lagrange-multiplier implementation fails the scalar kinetic gate. After integrating the multiplier term by parts, the velocity Hessian in `(Rdot, lambdadot)` is

```text
K = [[Z(R), ell^2], [ell^2, 0]]
det(K) = -ell^4
```

For any finite nonzero `ell`, one scalar kinetic eigenvalue is negative. The reduced background branch is retained only as a diagnostic and is not a physical solution of the full multiplier action.

## Structure

- `main.tex` — manuscript source.
- `graf.tex` — figure inclusion block, using `figures/` paths.
- `references.bib` — bibliography.
- `figures/` — generated figure assets.
- `code/locked_r_universe.py` — reduced-branch background diagnostic.
- `code/full_lambda_stability_gate.py` — full multiplier kinetic-gate calculation.
- `code/run_all.py` — reproduction plus static audit.
- `tables/` — sampled numerical outputs.

## Reproduce

```bash
python3 code/run_all.py
```
