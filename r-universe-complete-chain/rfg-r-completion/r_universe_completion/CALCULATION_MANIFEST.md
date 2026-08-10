# RFG-R Calculation Manifest

| Claim or calculation | Definition or file | Verification |
|---|---|---|
| Smooth RFG-R response | `scripts/rfg_regularized.py` | High-X and low-X tests |
| Exact potential reconstruction | `scripts/rfg_regularized.py` | Finite-difference ODE residual |
| Original cosmological branch recovery | `scripts/rfg_regularized.py` | Relative response and Q checks |
| Positive expanding background branch | `scripts/rfg_regularized.py` | Root residual and density closure |
| Tensor normalization and siren relation | `scripts/rfg_regularized.py` | `Q_T>0` on validation grid |
| Local GR matching function | `scripts/ppn_likelihood.py` | Solar Weyl ratio at one AU |
| PPN prediction in local domain | `scripts/ppn_likelihood.py` | `gamma=beta=1`, `alpha1=alpha2=0` |
| Cassini likelihood factor | `scripts/ppn_likelihood.py` | Direct Gaussian evaluation |
| Canonical-scalar sourced constraint reduction | `docs/canonical_scalar_reduction.md` | Exact Schur complement; GR limit `K=6`, `c_s^2=1` |
| Exact ADM-to-extended-EFT map | `scripts/extended_eft_mapping.py` | Independent implicit-background and `bar_m5` identity checks |
| Extended scalar-action audit | `scripts/extended_eft_scalar_stability.py` | Reduced action, lapse/shift determinant and scalar kinetic test |
| Full CMB/matter likelihood definition | `docs/likelihood_pipeline.md` | Blocked pending multi-fluid reduction after pure-gravity scalar degeneracy |
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
extended-EFT map           finite-difference checks and exact bar_m5 identity
extended scalar action     degenerate zeta-dot squared coefficient on 2,401-point (a,k) grid
```

The mapped coefficient `bar_m5=-M_Pl^2 Q_X/(3H0)` is nonzero on the reference
branch. The complete pure-gravity reduced extended-EFT action nevertheless has
a degenerate quadratic scalar kinetic term on the audited cosmological domain.
The physical multi-fluid reduction must determine whether matter lifts this
degeneracy or supplies an additional constraint before a CMB or matter
calculation. The public stock H-EFTCAMB coefficient interface also does not
expose `bar_m5`, so its compiled GR run remains a reference calculation only.
