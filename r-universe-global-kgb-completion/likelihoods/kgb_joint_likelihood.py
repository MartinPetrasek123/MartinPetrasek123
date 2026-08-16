"""Cobaya adapter for the executable covariant-KGB Planck plus late-time point."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

from cobaya.likelihood import Likelihood


class KGBJointLikelihood(Likelihood):
    """Slow, action-faithful likelihood with a cache per exact parameter point."""

    params = {
        "alpha": None,
        "omega_m0": None,
        "omega_r0": None,
        "ombh2": None,
        "logA": None,
        "ns": None,
        "tau": None,
        "A_planck": None,
    }

    source_root: str = ""
    output_root: str = ""
    python_executable: str = ""
    heftcamb_binary: str = ""
    heftcamb_template: str = ""
    clipy_source: str = ""
    planck_base: str = ""
    data_root: str = ""
    spline_points: int = 601
    late_time_integration_nodes: int = 16385

    def initialize(self) -> None:
        self._root = Path(self.source_root).resolve()
        self._evaluator = self._root / "scripts" / "evaluate_kgb_joint_point.py"
        self._output_root = Path(self.output_root).resolve()
        self._python = Path(getattr(self, "python_executable", sys.executable)).resolve()
        self._execution_python = Path(
            os.environ.get("KGB_EVALUATOR_PYTHON", str(self._python))
        ).resolve()
        for path, label in ((self._evaluator, "joint evaluator"), (self._python, "Python interpreter")):
            if not path.is_file():
                raise FileNotFoundError(f"{label} is unavailable: {path}")
        if not self._execution_python.is_file():
            raise FileNotFoundError(f"evaluator Python interpreter is unavailable: {self._execution_python}")
        provenance_files = {
            "joint_evaluator_sha256": self._evaluator,
            "generator_sha256": self._root / "scripts" / "generate_heftcamb_rph.py",
            "planck_wrapper_sha256": self._root / "scripts" / "evaluate_planck_2018_fixed.py",
            "late_time_evaluator_sha256": self._root / "scripts" / "evaluate_kgb_late_time.py",
            "heftcamb_binary_sha256": Path(self.heftcamb_binary).resolve(),
            "heftcamb_template_sha256": Path(self.heftcamb_template).resolve(),
        }
        missing = [label for label, path in provenance_files.items() if not path.is_file()]
        if missing:
            raise FileNotFoundError("missing evaluation-contract files: " + ", ".join(missing))
        # A cache entry is valid only for the exact executable calculation that
        # created it.  Paths pin the declared data release; hashes pin every
        # executable component that can change a numerical likelihood ordinate.
        self._execution_contract = {
            "format": 1,
            **{label: self._sha256(path) for label, path in provenance_files.items()},
            "clipy_source": str(Path(self.clipy_source).resolve()),
            "planck_base": str(Path(self.planck_base).resolve()),
            "data_root": str(Path(self.data_root).resolve()),
            "spline_points": int(self.spline_points),
            "late_time_integration_nodes": int(self.late_time_integration_nodes),
        }
        self._output_root.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    @staticmethod
    def _value(params: dict[str, Any], name: str) -> float:
        value = float(params[name])
        if not math.isfinite(value):
            raise ValueError(f"non-finite sampled parameter: {name}")
        return value

    def _point_key(self, inputs: dict[str, float]) -> str:
        payload = json.dumps(
            {"inputs": inputs, "execution_contract": self._execution_contract},
            sort_keys=True,
            separators=(",", ":"),
        ).encode("ascii")
        return hashlib.sha256(payload).hexdigest()[:24]

    def logp(self, **params_values: Any) -> float:
        alpha = self._value(params_values, "alpha")
        omega_m0 = self._value(params_values, "omega_m0")
        omega_r0 = self._value(params_values, "omega_r0")
        ombh2 = self._value(params_values, "ombh2")
        log_a_s = self._value(params_values, "logA")
        ns = self._value(params_values, "ns")
        tau = self._value(params_values, "tau")
        a_planck = self._value(params_values, "A_planck")
        scalar_amp = math.exp(log_a_s) * 1.0e-10
        inputs = {
            "alpha": alpha,
            "omega_m0": omega_m0,
            "omega_r0": omega_r0,
            "ombh2": ombh2,
            "scalar_amp": scalar_amp,
            "scalar_spectral_index": ns,
            "tau": tau,
            "a_planck": a_planck,
        }
        run_dir = self._output_root / self._point_key(inputs)
        command = [
            str(self._execution_python), str(self._evaluator),
            "--output-dir", str(run_dir),
            "--binary", str(self.heftcamb_binary),
            "--template", str(self.heftcamb_template),
            "--clipy-source", str(self.clipy_source),
            "--planck-base", str(self.planck_base),
            "--data-root", str(self.data_root),
            "--python", str(self._execution_python),
            "--alpha", format(alpha, ".17g"),
            "--omega-m0", format(omega_m0, ".17g"),
            "--omega-r0", format(omega_r0, ".17g"),
            "--ombh2", format(ombh2, ".17g"),
            "--scalar-amp", format(scalar_amp, ".17g"),
            "--scalar-spectral-index", format(ns, ".17g"),
            "--tau", format(tau, ".17g"),
            "--a-planck", format(a_planck, ".17g"),
            "--spline-points", str(self.spline_points),
            "--late-time-integration-nodes", str(self.late_time_integration_nodes),
            "--reuse-completed",
        ]
        if os.environ.get("KGB_COMPACT_ARTIFACTS") == "1":
            command.append("--compact-artifacts")
        result = subprocess.run(command, cwd=self._root, text=True, capture_output=True)
        if result.returncode:
            failure_log = run_dir / "cobaya_failure.log"
            failure_log.parent.mkdir(parents=True, exist_ok=True)
            failure_log.write_text(result.stdout + result.stderr, encoding="utf-8")
            return -math.inf
        summary_path = run_dir / "joint_point.json"
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            value = float(summary["minus_2_log_likelihood"]["total"])
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            return -math.inf
        return -0.5 * value if math.isfinite(value) else -math.inf
