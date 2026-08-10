#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/validate_completion.py
python3 scripts/canonical_scalar.py
python3 scripts/validate_extended_eft_mapping.py
python3 scripts/extended_eft_mapping.py
python3 scripts/validate_extended_eft_scalar_stability.py
python3 scripts/validate_multifluid_reduction.py
python3 scripts/generate_outputs.py
python3 scripts/ppn_likelihood.py
echo "RFG-R action, exact finite multi-fluid reduction, and hierarchy interface checks passed; no RFG-R data likelihood is evaluated by this command."
