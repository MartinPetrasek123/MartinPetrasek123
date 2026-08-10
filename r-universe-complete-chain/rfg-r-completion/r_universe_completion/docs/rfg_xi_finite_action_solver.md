# RFG-RXi Finite-Action DAE Solver Kernel

`scripts/rfg_xi_finite_action_solver.py` is the first executable time-domain
solver generated directly from the completed RFG-RXi quadratic action. It is
not a CMB-spectrum or likelihood solver. Its purpose is narrower and
testable: demonstrate that the action-level material DAE can be integrated
from `a=1e-8` to `a=1` without replacing its metric constraints by Einstein
equations or a quasi-static effective coupling.

For the reduced variables `y=(zeta, Delta_i)`, the action gives `K_0i=0` and a
symmetric mixing matrix. The curvature row is therefore the algebraic
constraint

```text
mu_zeta zeta + [M_0i-Bdot_0i] Delta_i = 0.
```

The material rows are integrated in `N=ln(a)` as

```text
K_ij [H^2 Delta_j,NN + Hdot Delta_j,N]
 + Kdot_ij H Delta_j,N
 + [Bdot_iA-M_iA] y_A = 0.
```

Every `K`, `Kdot`, `Bdot`, and `M` entry is evaluated from the analytic ADM
action derivatives. The code obtains `zeta` by solving the algebraic row at
each right-hand-side evaluation and rejects a vanishing normalized
`mu_zeta`.

The current kernel contains the four material monopole/dipole variables that
belong to the exact finite action. Photon polarization, photon higher moments,
collisionless neutrino higher moments, recombination, line-of-sight sources,
lensing, and primordial adiabatic normalization are deliberately not yet
added; consequently its CSV trajectory is not labelled a transfer function or
a spectrum.

Run:

```bash
cd r_universe_completion
python3 -m pip install scipy==1.13.1
python3 scripts/rfg_xi_finite_action_solver.py
python3 scripts/validate_rfg_xi_finite_action_solver.py
```

The generated table is `generated/tables/rfg_xi_finite_action_solver.csv`.
The regression test checks finite RFG-RXi trajectories for Xi=1 and Xi=2,
constraint residual below `1e-16`, exact linear amplitude scaling, and the
Einstein-Hilbert action limit `Q=1`, `Q_X=m5_bar=0`.
