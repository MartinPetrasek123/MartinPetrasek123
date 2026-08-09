#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/validate_completion.py
python3 scripts/generate_outputs.py
python3 scripts/ppn_likelihood.py
echo "All RFG-R completion checks passed."
