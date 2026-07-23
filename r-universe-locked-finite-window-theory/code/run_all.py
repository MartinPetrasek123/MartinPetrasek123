#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CODEX_PY = Path('/Users/mpetr/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3')
PY = os.environ.get('RUN_ALL_PY') or (str(DEFAULT_CODEX_PY) if DEFAULT_CODEX_PY.exists() else sys.executable)


def static_audit() -> None:
    tex = (ROOT / 'main.tex').read_text() + (ROOT / 'graf.tex').read_text()
    labels = re.findall(r'\\label\{([^}]+)\}', tex)
    refs = re.findall(r'\\(?:ref|eqref)\{([^}]+)\}', tex)
    imgs = re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', (ROOT / 'graf.tex').read_text())
    forbidden = ['LCDM', 'Pantheon', 'DESI', 'benchmark', 'against', 'relative to', 'standard cosmology', 'kgb_reconstruction', 'Enter Caption', 'fig:placeholder', 'Fig. ??']
    out = {
        'missing_refs': sorted(set(refs) - set(labels)),
        'duplicate_labels': [k for k, v in Counter(labels).items() if v > 1],
        'includegraphics': imgs,
        'missing_graphics_in_root': [img for img in imgs if not (ROOT / img).exists()],
        'forbidden_phrase_hits': {s: tex.count(s) for s in forbidden},
    }
    (ROOT / 'code' / 'latex_static_audit.json').write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(json.dumps(out, indent=2, ensure_ascii=False))


def main() -> None:
    start = time.time()
    subprocess.run([PY, str(ROOT / 'code' / 'locked_r_universe.py')], cwd=ROOT, check=True)
    subprocess.run([PY, str(ROOT / 'code' / 'full_lambda_stability_gate.py')], cwd=ROOT, check=True)
    static_audit()
    print(f'[run_all] locked R-Universe reproduction completed in {time.time() - start:.2f} s')


if __name__ == '__main__':
    main()
