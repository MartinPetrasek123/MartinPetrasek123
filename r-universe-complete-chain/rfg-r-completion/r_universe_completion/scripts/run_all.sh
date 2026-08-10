#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/validate_completion.py
python3 scripts/canonical_scalar.py
python3 scripts/validate_extended_eft_mapping.py
python3 scripts/extended_eft_mapping.py
python3 scripts/validate_extended_eft_scalar_stability.py
python3 scripts/validate_multifluid_reduction.py
python3 scripts/rfg_dae_closure.py
python3 scripts/validate_rfg_dae_closure.py
python3 scripts/rfg_xi_completion.py
python3 scripts/validate_rfg_xi_completion.py
python3 scripts/validate_rfg_xi_metric_equations.py
python3 scripts/rfg_xi_observables.py
python3 scripts/validate_rfg_xi_observables.py
python3 scripts/generate_outputs.py
python3 scripts/ppn_likelihood.py
python3 scripts/validate_manuscript_links.py --manuscript ../manuscript/main.tex --package-root .
echo "RFG-R and RFG-RXi action checks passed; neither model has an executed CMB/matter data likelihood."
