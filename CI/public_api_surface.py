#!/usr/bin/env python3
"""Measure the ACTS public API surface from Doxygen XML.

The "public API surface" is the set of *documented*, non-``detail`` /
non-``Experimental`` entities in the ``Acts`` namespace, matching the policy in
docs/pages/versioning.md. Doxygen parses syntactically, so this needs no
dependency headers or compile flags -- just the ``doxygen`` binary.

Typical use (runs Doxygen itself, then reports):

    CI/public_api_surface.py --run --json surface.json --markdown surface.md

Or parse an existing XML tree:

    CI/public_api_surface.py --xml <dir>/xml --markdown -

Writes a Markdown summary to ``$GITHUB_STEP_SUMMARY`` when that env var is set.
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

DOXYFILE = "CI/public_api/Doxyfile"
INPUT_ROOT = "Core/include/Acts"

# Doxygen memberdef 'kind' -> report bucket for namespace-scope members
NS_MEMBER_BUCKET = {
    "function": "free_functions",
    "typedef": "aliases",
    "variable": "variables",
    "enum": "enums",
    "concept": "concepts",
}


def is_internal(name: str) -> bool:
    return bool(name) and ("detail" in name or "Experimental" in name)


def module_of(location: str | None) -> str:
    """Top-level module from a location like 'Core/include/Acts/Surfaces/X.hpp'."""
    if not location:
        return "?"
    parts = Path(location).as_posix().split("/")
    if "Acts" in parts:
        i = parts.index("Acts")
        if i + 1 < len(parts) - 1:  # something between Acts/ and the filename
            return parts[i + 1]
    return "(root)"


def run_doxygen(repo: Path, out: Path, input_dir: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, DOXY_OUT=str(out), DOXY_INPUT=str(input_dir))
    subprocess.run(["doxygen", DOXYFILE], cwd=repo, env=env, check=True)


def _norm_type(el) -> str:
    """Flatten a Doxygen <type>/<param><type> element to normalized text."""
    if el is None:
        return ""
    return " ".join("".join(el.itertext()).split())


def callable_forms(md, qualname: str) -> dict[str, str]:
    """Expand a function memberdef into its source-callable signatures.

    A parameter with a default value makes shorter calls valid too, so
    ``f(A, B = d)`` yields both ``f(A)`` and ``f(A, B)``. Comparing these
    expanded forms across two revisions makes *adding a defaulted argument*
    non-breaking while *adding a non-defaulted argument* (or removing/retyping
    one) shows up as a removed form. Maps form-key -> return type.
    """
    params = md.findall("param")
    types = [_norm_type(p.find("type")) for p in params]
    has_default = [p.find("defval") is not None for p in params]
    required = len(types)
    for i, d in enumerate(has_default):
        if d:
            required = i
            break
    ret = _norm_type(md.find("type"))
    const = " const" if md.get("const") == "yes" else ""
    forms = {}
    for m in range(required, len(types) + 1):
        forms[f"{qualname}({', '.join(types[:m])}){const}"] = ret
    return forms


def parse_xml(xml_dir: Path) -> dict:
    counts: Counter[str] = Counter()
    per_module: dict[str, Counter] = defaultdict(Counter)
    type_names: set[str] = set()
    ns_member_names: set[str] = set()
    symbols: set[str] = set()  # stable per-symbol keys for diffing (non-function)
    callables: dict[str, str] = {}  # source-callable signature -> return type
    methods = 0

    for f in glob.glob(str(xml_dir / "*.xml")):
        base = os.path.basename(f)
        if base in ("index.xml", "Doxyfile.xml"):
            continue
        try:
            root = ET.parse(f).getroot()
        except ET.ParseError:
            continue
        for cd in root.findall("compounddef"):
            kind = cd.get("kind")
            name = cd.findtext("compoundname") or ""
            if not name.startswith("Acts") or is_internal(name):
                continue
            loc = cd.find("location")
            locfile = loc.get("file") if loc is not None else None

            if kind in ("class", "struct", "union"):
                bare = name.split("<")[0]
                if bare not in type_names:
                    type_names.add(bare)
                    counts["types"] += 1
                    per_module[module_of(locfile)]["types"] += 1
                    symbols.add(f"type {bare}")
                # public methods on this documented type
                for md in cd.iter("memberdef"):
                    if md.get("kind") == "function" and md.get("prot") == "public":
                        methods += 1
                        callables.update(
                            callable_forms(md, f"{bare}::{md.findtext('name') or ''}"))

            elif kind == "concept":
                # Doxygen >= 1.9.7 emits C++20 concepts as their own compound.
                if name not in ns_member_names:
                    ns_member_names.add(name)
                    counts["concepts"] += 1
                    per_module[module_of(locfile)]["concepts"] += 1
                    symbols.add(f"concept {name}")

            elif kind == "namespace":
                for md in cd.findall("sectiondef/memberdef"):
                    mk = md.get("kind")
                    bucket = NS_MEMBER_BUCKET.get(mk)
                    if not bucket:
                        continue
                    mname = md.findtext("name") or ""
                    full = f"{name}::{mname}"
                    if full in ns_member_names:
                        continue
                    ns_member_names.add(full)
                    counts[bucket] += 1
                    symbols.add(f"{bucket} {full}")
                    if mk == "function":
                        callables.update(callable_forms(md, full))
                    mloc = md.find("location")
                    per_module[module_of(mloc.get("file") if mloc is not None else None)][bucket] += 1

    total = sum(counts.values())
    return {
        "counts": dict(counts),
        "total": total,
        "public_methods": methods,
        "per_module": {m: dict(c) for m, c in per_module.items()},
        "symbols": sorted(symbols),
        "callables": callables,
    }


BUCKET_ORDER = ["types", "free_functions", "aliases", "variables", "enums", "concepts"]


def render_markdown(data: dict, doxy_version: str | None) -> str:
    c = data["counts"]
    lines = ["## ACTS public API surface", ""]
    lines.append(f"**{data['total']}** documented public names in `Acts::` "
                 f"(excluding `detail` / `Experimental`), plus "
                 f"**{data['public_methods']}** public methods on documented types.")
    if doxy_version:
        lines.append(f"\n_Source: `Core/include/Acts` via Doxygen {doxy_version}._")
    lines += ["", "| category | count |", "|---|---:|"]
    for k in BUCKET_ORDER:
        if k in c:
            lines.append(f"| {k.replace('_', ' ')} | {c[k]} |")
    lines.append(f"| **total** | **{data['total']}** |")

    lines += ["", "<details><summary>By module</summary>", "",
              "| module | total | " + " | ".join(BUCKET_ORDER) + " |",
              "|---|---:|" + "|".join(["---:"] * len(BUCKET_ORDER)) + "|"]
    mods = sorted(data["per_module"].items(), key=lambda kv: -sum(kv[1].values()))
    for m, d in mods:
        row = [m, str(sum(d.values()))] + [str(d.get(k, 0)) for k in BUCKET_ORDER]
        lines.append("| " + " | ".join(row) + " |")
    lines += ["", "</details>"]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", default=".", help="repository root")
    ap.add_argument("--run", action="store_true", help="run doxygen (else use --xml)")
    ap.add_argument("--input", help="header root to measure "
                    "(default <repo>/Core/include/Acts); use to point at another checkout")
    ap.add_argument("--xml", help="existing Doxygen XML dir to parse")
    ap.add_argument("--json", help="write report JSON here")
    ap.add_argument("--markdown", help="write Markdown here ('-' for stdout)")
    args = ap.parse_args()

    repo = Path(args.repo).resolve()

    doxy_version = None
    try:
        doxy_version = subprocess.run(["doxygen", "--version"],
                                      capture_output=True, text=True).stdout.split()[0]
    except (FileNotFoundError, IndexError):
        pass

    tmp = None
    if args.run:
        input_dir = Path(args.input).resolve() if args.input else repo / INPUT_ROOT
        tmp = Path(tempfile.mkdtemp(prefix="acts-api-surface-"))
        run_doxygen(repo, tmp, input_dir)
        xml_dir = tmp / "xml"
    elif args.xml:
        xml_dir = Path(args.xml)
    else:
        print("error: pass --run or --xml", file=sys.stderr)
        return 2

    if not xml_dir.is_dir():
        print(f"error: no XML at {xml_dir}", file=sys.stderr)
        return 2

    data = parse_xml(xml_dir)
    data["doxygen_version"] = doxy_version

    md = render_markdown(data, doxy_version)

    if args.json:
        import json
        Path(args.json).write_text(json.dumps(data, indent=2) + "\n")
    if args.markdown == "-" or args.markdown is None:
        print(md)
    elif args.markdown:
        Path(args.markdown).write_text(md)

    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as fh:
            fh.write(md)

    return 0


if __name__ == "__main__":
    sys.exit(main())
