#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
python3 scripts/validate_ru_kgb.py
python3 scripts/stability_scan.py
python3 scripts/matter_qs.py
python3 scripts/ppn_screening.py
python3 scripts/cmb_prerecombination.py
python3 scripts/make_figures.py
echo "Global R-alpha KGB reconstruction and stability gates passed."
