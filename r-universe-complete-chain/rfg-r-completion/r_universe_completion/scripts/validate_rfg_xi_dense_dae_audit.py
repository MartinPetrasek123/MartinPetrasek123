#!/usr/bin/env python3
"""Validate the recorded extended RFG-RXi DAE audit without rerunning it."""

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "generated" / "tables" / "rfg_xi_dense_dae_audit.csv"


def main() -> None:
    with TABLE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 386
    assert {float(row["xi_completion"]) for row in rows} == {1.0, 2.0}
    for xi_completion in (1.0, 2.0):
        branch = [row for row in rows if float(row["xi_completion"]) == xi_completion]
        assert len(branch) == 193
        assert math.isclose(float(branch[0]["a"]), 1.0e-8, rel_tol=0.0, abs_tol=1.0e-22)
        assert math.isclose(float(branch[-1]["a"]), 1.0, rel_tol=0.0, abs_tol=1.0e-14)
        assert all(row["root_k_over_H0"] == "" for row in branch)
        assert min(float(row["min_mu_zeta_relative_on_k_grid"]) for row in branch) > 1.0e-8
        assert min(float(row["min_Q_tensor_on_k_grid"]) for row in branch) > 0.0
        assert min(float(row["min_constraint_determinant_on_k_grid"]) for row in branch) > 0.0
    print("RFG-RXi dense DAE audit regression checks passed")


if __name__ == "__main__":
    main()
