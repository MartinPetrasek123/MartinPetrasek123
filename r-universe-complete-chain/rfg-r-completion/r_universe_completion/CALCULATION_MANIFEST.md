# RFG-R Calculation Manifest

| Claim or calculation | Definition or file | Verification |
|---|---|---|
| Smooth RFG-R response | `scripts/rfg_regularized.py` | High-X and low-X tests |
| Exact potential reconstruction | `scripts/rfg_regularized.py` | Defining reconstruction identity; finite-difference residual used only as an independent regression check |
| Original cosmological branch recovery | `scripts/rfg_regularized.py` | Relative response and Q checks |
| Positive expanding background branch | `scripts/rfg_regularized.py` | Root residual and density closure |
| Tensor normalization and siren relation | `scripts/rfg_regularized.py` | `Q_T>0` on validation grid |
| Local GR matching function | `scripts/ppn_likelihood.py` | Solar Weyl ratio at one AU |
| PPN prediction in local domain | `scripts/ppn_likelihood.py` | `gamma=beta=1`, `alpha1=alpha2=0` |
| Cassini likelihood factor | `scripts/ppn_likelihood.py` | Direct Gaussian evaluation |
| Canonical-scalar sourced constraint reduction | `docs/canonical_scalar_reduction.md` | Exact Schur complement; GR limit `K=6`, `c_s^2=1` |
| Exact ADM-to-extended-EFT map | `scripts/extended_eft_mapping.py` | Independent implicit-background and `bar_m5` identity checks |
| Extended scalar-action audit | `scripts/extended_eft_scalar_stability.py` | Reduced action, lapse/shift determinant and scalar kinetic test |
| Photon--baryon--CDM--neutrino reduction | `docs/photon_baryon_cdm_neutrino_reduction.md` | Exact Sorkin--Schutz Schur complement, infinite photon/neutrino hierarchy, GR rational-rank and RFG-R inertia audits |
| Multi-fluid DAE closure | `scripts/rfg_dae_closure.py` | Exact `Delta_i` transformation, analytic Schur-complement time derivative, and singular curvature-constraint surface |
| RFG-RXi completion | `scripts/rfg_xi_completion.py`, `scripts/validate_rfg_xi_metric_equations.py`, `scripts/rfg_xi_dense_dae_audit.py` | Directly derived Xi(R3+sigma^2) increments in the lapse, shift, trace, and traceless equations; unchanged FLRW background; no-root 49-by-81 audit; and extended no-root 193-by-601 audit at unfitted Xi=1 and Xi=2 |
| RFG-RXi direct data factors | `scripts/rfg_xi_observables.py` | Exact background/tensor observables plus the locally matched Cassini likelihood for Xi=1 and Xi=2; explicitly not a CMB likelihood |
| GR CAMB spectrum/transfer reference | `scripts/gr_reference_camb.py` | Pinned CAMB 2.0.1 regression for lensed CMB, lensing potential and linear matter transfer; never an RFG-R spectrum |
| GR Planck low-ell+lensing interface | `scripts/gr_planck_lowell_lensing.py` | Executed official low-T, low-E and lensing likelihood regression at the pinned GR point; high-ell nuisance vector not guessed |
| Full CMB/matter likelihood definition | `docs/likelihood_pipeline.md` | Exact action and hierarchy interface specified; no RFG-R or RFG-RXi spectra or data likelihood has been evaluated |
| Standalone article | `paper/R_Universe_RFG_R_Completion.pdf` | Rendered PDF and text extraction |

## Reference Validation

Run:

```bash
bash scripts/run_all.sh
```

Reference results:

```text
max high-X response error  <= 1e-12 (machine precision in practice)
max high-X Q error         <= 1e-12 (machine precision in practice)
max potential ODE residual <= 2e-5
max relative branch error  <= 1e-12
max density closure error  <= 1e-11
min Q_T                    > 0
Solar W at one AU          > 1e9
canonical scalar GR limit  K=6, c_s^2=1
extended-EFT map           analytic action derivatives; finite-difference regression and exact bar_m5 identity
extended scalar action     degenerate zeta-dot squared coefficient on 2,401-point (a,k) grid
multi-fluid core           exact GR rank(K)=4; RFG-R reference inertia (4,0,1) on 425 points
multi-fluid DAE            zeta is algebraically constrained, but mu_zeta=0 has 24/49 root-bearing scale factors; a=1 root k/H0=2.51545672221
RFG-RXi DAE               unfitted Xi=1 and Xi=2 have no sampled mu_zeta root on separate 49x81 (a,k) grids; no cosmological data likelihood has been run
RFG-RXi direct factors    Cassini -2 ln L=0.8336483931947068; c_T-1=0; dL_GW/dL_EM(z~1)=0.9565968846 (Xi=1), 0.9717969466 (Xi=2)
GR CAMB reference          CAMB 2.0.1, sigma8(z=0)=0.8110325278646
GR Planck low-ell+lensing  -2 ln L=428.3415086187 at the fixed GR reference point
```

The mapped coefficient `bar_m5=-M_Pl^2 Q_X/(3H0)` is nonzero on the reference
branch. The complete pure-gravity reduced extended-EFT action nevertheless has
a degenerate quadratic scalar kinetic term on the audited cosmological domain.
The exact finite multi-fluid reduction supplies one null constraint and four
positive monopole/dipole material directions; photons and neutrinos retain
their untruncated kinetic hierarchy. The exact DAE closure identifies that
null direction with the curvature constraint and finds its coefficient crosses
zero on the reference branch. This blocks a globally regular RFG-R spectrum
or likelihood for that branch. RFG-RXi is a separately stated completion with
the background-null operator Xi(R3+sigma^2); it passes the present finite-grid
DAE root audit, but has no CMB/matter inference. The audited public stock
H-EFTCAMB coefficient interface also does not expose `bar_m5`, so any stock
GR run remains a reference calculation only.
