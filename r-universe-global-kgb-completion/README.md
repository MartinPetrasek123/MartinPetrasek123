# R-Universe global R-alpha KGB completion

This package turns the R-alpha background already fitted in the accompanying
late-time R-Universe manuscript into one explicit, globally defined,
luminal KGB action.  It replaces the missing `kgb_reconstruction.py` named by
that manuscript with reproducible code and closes the scalar no-ghost and
no-gradient conditions analytically and numerically.

Run:

```sh
bash scripts/run_all.sh
```

Outputs:

- `generated/ru_kgb_trajectory.csv` - action and EFT functions;
- `generated/validation.json` - all numerical gates;
- `generated/stability_scan.csv` and `.json` - physical stability scan around the fit;
- `generated/matter_qs_prediction.csv` - derived high-k linear matter response and growth;
- `generated/ppn_screening.json` - solar-system cubic-KGB screening gate;
- `generated/cmb_prerecombination_gate.json` - sound-horizon and pre-recombination EFT inputs;
- `generated/ru_kgb_stability.pdf` and `.png` - stability figure.

The standalone technical paper source is
`paper/R_Universe_Global_KGB_Completion.tex`. It references this exact source
tree and can be included as the covariant-completion companion to the
late-time R-Universe manuscript.

The construction preserves the R-alpha likelihood background exactly for all
observed redshifts.  It includes a derived subhorizon matter prediction and a
decoupling-limit PPN screening gate. A CMB angular-spectrum likelihood remains
a separate executable Boltzmann calculation rather than a fluid replacement.
