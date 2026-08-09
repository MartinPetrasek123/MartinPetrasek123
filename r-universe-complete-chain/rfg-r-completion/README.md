# RFG-R complete source and reproducibility package

This directory contains the complete LaTeX source of the RCD/RFG/RFG-R paper
and the executable RFG-R calculation package it cites.

The public canonical location is
<https://github.com/MartinPetrasek123/MartinPetrasek123/tree/main/r-universe-complete-chain/rfg-r-completion>.

## Read in this order

1. `main.tex` -- full manuscript. It contains the primary action, the exact
   background reconstruction, the controlled recovery limit, the local GR
   matching rule, the PPN prediction, the tensor result, and the exact CMB
   likelihood contract.
2. `r_universe_completion/docs/completion_derivation.md` -- complete RFG-R
   action-level derivation and limits.
3. `r_universe_completion/CALCULATION_MANIFEST.md` -- claim-to-file mapping.
4. `r_universe_completion/docs/likelihood_pipeline.md` -- required full
   matter/CMB/PPN data-likelihood calculation, including rejection rules.
5. `r_universe_completion/scripts/run_all.sh` -- deterministic local
   validation and output generation.

## Reproduce the internal RFG-R checks

From this directory:

```bash
bash r_universe_completion/scripts/run_all.sh
```

The command validates regularity, recovery of the original cosmological
branch, potential reconstruction, background closure, positive tensor kinetic
normalization, local GR matching, and the Cassini likelihood factor. It also
regenerates the committed tables and figures.

## Scientific status

RFG-R is a defined low-energy multiscale EFT, not a declaration that it has
already replaced LambdaCDM. Its local PPN prediction follows from the stated
exact-GR matching domain. Its matter/CMB calculation is a complete
implementation and likelihood specification, but no compiled RFG-R
Einstein--Boltzmann module or joint Planck/BAO/SN/RSD/GW posterior is included
yet. The public package and manuscript state this boundary explicitly.

## Stable public links

- [Manuscript source](main.tex)
- [Action derivation](r_universe_completion/docs/completion_derivation.md)
- [Calculation manifest](r_universe_completion/CALCULATION_MANIFEST.md)
- [Matter/CMB/PPN likelihood protocol](r_universe_completion/docs/likelihood_pipeline.md)
- [Core model code](r_universe_completion/scripts/rfg_regularized.py)
- [Validation code](r_universe_completion/scripts/validate_completion.py)
- [PPN likelihood code](r_universe_completion/scripts/ppn_likelihood.py)
- [Generated tables and figures](r_universe_completion/generated)
- [Standalone RFG-R paper](r_universe_completion/paper/R_Universe_RFG_R_Completion.pdf)

For a citable release, create an immutable GitHub release and archive that
release with Zenodo before journal submission. Do not cite a DOI until Zenodo
has issued it.
