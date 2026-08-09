#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/validate_completion.py
python3 scripts/canonical_scalar.py
python3 scripts/generate_outputs.py
python3 scripts/ppn_likelihood.py
echo "Conditional RFG-R branch checks passed; no global-completion claim is made."
