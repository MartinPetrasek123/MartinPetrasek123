#!/usr/bin/env python3
"""Independent regression checks for directly testable RFG-RXi observables."""

import math

from rfg_xi_completion import XI_AUDIT_VALUES
from rfg_xi_observables import xi_observable_rows


def main() -> None:
    rows = xi_observable_rows()
    expected_count = len(XI_AUDIT_VALUES) * 301
    assert len(rows) == expected_count
    for xi_completion in XI_AUDIT_VALUES:
        branch = [row for row in rows if row["xi_completion"] == xi_completion]
        today = branch[-1]
        assert today["a"] == 1.0
        assert math.isclose(today["dL_GW_over_dL_EM"], 1.0, rel_tol=0.0, abs_tol=1.0e-14)
        assert today["Q_tensor"] > 0.0
        assert all(row["c_T"] == 1.0 for row in branch)
        assert all(row["gamma_PPN"] == 1.0 and row["beta_PPN"] == 1.0 for row in branch)
        assert all(row["alpha1_PPN"] == 0.0 and row["alpha2_PPN"] == 0.0 for row in branch)
        assert all(row["Q_tensor"] > 0.0 for row in branch)
    first = [row for row in rows if row["xi_completion"] == XI_AUDIT_VALUES[0]]
    second = [row for row in rows if row["xi_completion"] == XI_AUDIT_VALUES[1]]
    assert len(first) == len(second)
    for left, right in zip(first, second):
        for key in ("E", "H_km_s_Mpc", "Omega_m", "Omega_r", "Omega_R"):
            assert left[key] == right[key]
    print("RFG-RXi observable regression checks passed")
    print("background is exactly Xi-independent; tensor and local PPN checks are positive")


if __name__ == "__main__":
    main()
