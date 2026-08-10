#!/usr/bin/env python3
"""Check that the public LaTeX manuscript is internally linked and reproducible.

This is deliberately a source-level gate: it is independent of a TeX binary,
so a missing local TeX installation cannot hide a dangling scientific link.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen


LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
REF_RE = re.compile(r"\\(?:auto|eq|page|v|V|c|C)?ref\{([^}]+)\}")
CITE_RE = re.compile(r"\\cite[a-zA-Z*]*\{([^}]+)\}")
BIB_RE = re.compile(r"\\bibitem(?:\[[^]]*\])?\{([^}]+)\}")
GRAPHIC_RE = re.compile(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}")
REPO_RE = re.compile(r"\\RFGRepo(?:Named|Path)\{([^}]+)\}")
REPO_URL_RE = re.compile(r"\\newcommand\{\\RFGRepoURL\}\{([^}]+)\}")
ENVIRONMENT_RE = re.compile(r"\\(begin|end)\{([A-Za-z*]+)\}")
MATH_COMMAND_RE = re.compile(
    r"\\(?:alpha|beta|gamma|delta|epsilon|eta|theta|zeta|Xi|Pi|rho|"
    r"sigma|mu|nu|Omega|dot|ddot|frac|sqrt|sum|int|partial|mathrm|rm|"
    r"mathcal|mathbf|operatorname)\b"
)
MATH_ENVIRONMENTS = {
    "align",
    "align*",
    "aligned",
    "alignedat",
    "array",
    "bmatrix",
    "cases",
    "displaymath",
    "equation",
    "equation*",
    "gather",
    "gather*",
    "math",
    "matrix",
    "multline",
    "multline*",
    "pmatrix",
    "smallmatrix",
    "split",
    "Vmatrix",
    "vmatrix",
}


def _tex_path(value: str) -> Path:
    """Convert the restricted path syntax used inside the manuscript."""
    decoded = value.replace(r"\_", "_").replace(r"\-", "-")
    path = Path(decoded)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe manuscript path: {value}")
    return path


def _split_keys(matches: list[str]) -> set[str]:
    return {
        key.strip()
        for value in matches
        for key in value.split(",")
        if key.strip()
    }


def _local_figure_path(manuscript: Path, value: str) -> Path | None:
    relative = _tex_path(value)
    candidates = [manuscript.parent / relative]
    if not relative.suffix:
        candidates.extend(
            manuscript.parent / relative.with_suffix(suffix)
            for suffix in (".pdf", ".png", ".jpg", ".jpeg")
        )
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def _remote_check(base_url: str, paths: list[Path]) -> list[str]:
    failures: list[str] = []
    prefix = "r-universe-complete-chain/rfg-r-completion/r_universe_completion/"
    for path in paths:
        url = f"{base_url}/blob/main/{quote(prefix + path.as_posix(), safe='/')}"
        try:
            request = Request(url, headers={"User-Agent": "r-universe-manuscript-check/1.0"})
            with urlopen(request, timeout=20) as response:
                if response.status != 200:
                    failures.append(f"HTTP {response.status}: {url}")
        except Exception as exc:  # noqa: BLE001 - test must report transport failures
            failures.append(f"{url}: {exc}")
    return failures


def _environment_errors(source: str) -> list[str]:
    """Return nesting errors for LaTeX environments without requiring TeX."""
    stack: list[str] = []
    errors: list[str] = []
    for action, environment in ENVIRONMENT_RE.findall(source):
        if action == "begin":
            stack.append(environment)
        elif not stack:
            errors.append(f"end of unopened environment: {environment}")
        elif stack[-1] != environment:
            errors.append(
                f"environment mismatch: expected end of {stack[-1]}, got {environment}"
            )
        else:
            stack.pop()
    errors.extend(f"unclosed environment: {environment}" for environment in reversed(stack))
    return errors


def _math_mode_errors(source: str) -> list[str]:
    """Find mathematical tokens that LaTeX would reject in ordinary prose."""
    document_marker = r"\begin{document}"
    start = source.find(document_marker)
    body = source[start + len(document_marker) :] if start >= 0 else source
    errors: list[str] = []
    math_depth = 0
    dollar_mode: str | None = None
    i = 0
    while i < len(body):
        if body[i] == "%" and (i == 0 or body[i - 1] != "\\"):
            newline = body.find("\n", i)
            i = len(body) if newline < 0 else newline + 1
            continue
        begin_end = ENVIRONMENT_RE.match(body, i)
        if begin_end:
            action, environment = begin_end.groups()
            if environment in MATH_ENVIRONMENTS:
                math_depth += 1 if action == "begin" else -1
            i = begin_end.end()
            continue
        if body.startswith(r"\(", i) or body.startswith(r"\[", i):
            math_depth += 1
            i += 2
            continue
        if body.startswith(r"\)", i) or body.startswith(r"\]", i):
            math_depth = max(0, math_depth - 1)
            i += 2
            continue
        if body[i] == "$" and (i == 0 or body[i - 1] != "\\"):
            marker = "$$" if body.startswith("$$", i) else "$"
            if dollar_mode is None:
                dollar_mode = marker
                math_depth += 1
            elif dollar_mode == marker:
                dollar_mode = None
                math_depth = max(0, math_depth - 1)
            i += len(marker)
            continue
        if math_depth == 0:
            if body[i] in "_^" and (i == 0 or body[i - 1] != "\\"):
                line = body.count("\n", 0, i) + 1
                errors.append(f"raw '{body[i]}' outside math mode near document line {line}")
            command = MATH_COMMAND_RE.match(body, i)
            if command:
                line = body.count("\n", 0, i) + 1
                errors.append(f"math command {command.group(0)} outside math mode near document line {line}")
                i = command.end()
                continue
        i += 1
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", type=Path, required=True)
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--check-remote", action="store_true")
    args = parser.parse_args()

    manuscript = args.manuscript.resolve()
    package_root = args.package_root.resolve()
    source = manuscript.read_text(encoding="utf-8")

    labels = LABEL_RE.findall(source)
    duplicate_labels = sorted(key for key, count in Counter(labels).items() if count > 1)
    referenced = _split_keys(REF_RE.findall(source))
    missing_refs = sorted(referenced - set(labels))

    citations = _split_keys(CITE_RE.findall(source))
    bibliography = set(BIB_RE.findall(source))
    missing_citations = sorted(citations - bibliography)

    figures = GRAPHIC_RE.findall(source)
    missing_figures = [figure for figure in figures if _local_figure_path(manuscript, figure) is None]
    environment_errors = _environment_errors(source)
    math_errors = _math_mode_errors(source)

    repo_paths = [_tex_path(value) for value in REPO_RE.findall(source)]
    missing_repo_paths = [path.as_posix() for path in repo_paths if not (package_root / path).is_file()]

    failures: list[str] = []
    if duplicate_labels:
        failures.append("duplicate labels: " + ", ".join(duplicate_labels))
    if missing_refs:
        failures.append("unresolved references: " + ", ".join(missing_refs))
    if missing_citations:
        failures.append("unresolved citations: " + ", ".join(missing_citations))
    if missing_figures:
        failures.append("missing figures: " + ", ".join(missing_figures))
    if environment_errors:
        failures.append("environment structure: " + "; ".join(environment_errors))
    if math_errors:
        failures.append("math-mode structure: " + "; ".join(math_errors))
    if missing_repo_paths:
        failures.append("missing linked calculations: " + ", ".join(missing_repo_paths))

    if args.check_remote:
        base_match = REPO_URL_RE.search(source)
        if base_match is None:
            failures.append("RFGRepoURL command not found")
        else:
            failures.extend(_remote_check(base_match.group(1), repo_paths))

    if failures:
        print("manuscript link validation failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print(
        "manuscript link validation passed: "
        f"{len(labels)} labels, {len(referenced)} references, "
        f"{len(citations)} citations, {len(figures)} figures, "
        f"{len(repo_paths)} linked calculation paths"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
