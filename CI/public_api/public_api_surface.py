#!/usr/bin/env python3
"""Measure the ACTS public API surface from Doxygen XML.

The "public API surface" is the set of *documented*, non-``detail`` /
non-``Experimental`` entities in the ``Acts`` namespace, matching the policy in
docs/pages/versioning.md. Doxygen parses syntactically, so this needs no
dependency headers or compile flags -- just the ``doxygen`` binary.

Typical use (runs Doxygen itself, then reports):

    CI/public_api/public_api_surface.py --run --json surface.json --markdown surface.md

Or parse an existing XML tree:

    CI/public_api/public_api_surface.py --xml <dir>/xml --markdown -

Writes a Markdown summary to ``$GITHUB_STEP_SUMMARY`` when that env var is set.
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import tempfile
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

# Absolute so it resolves regardless of the caller's working directory.
DOXYFILE = str(Path(__file__).resolve().parent / "Doxyfile")

# Plugins deliberately left out of the public-API surface metric. (The
# top-level Detray/ and Traccc/ build-integration folders live outside Plugins/
# and are never scanned; the Detray *plugin* under Plugins/ is in scope.)
EXCLUDED_PLUGINS: set[str] = set()


def standard_roots(under: Path) -> list[Path]:
    """Header roots whose public API we track, resolved under `under`.

    Core, Fatras and Alignment, plus every plugin under Plugins/. The Examples
    tree and the top-level Detray/Traccc integration folders are not included.
    """
    roots: list[Path] = []
    for rel in ("Core/include", "Fatras/include", "Alignment/include"):
        p = under / rel
        if p.is_dir():
            roots.append(p)
    plugins = under / "Plugins"
    if plugins.is_dir():
        for pdir in sorted(plugins.iterdir()):
            if (
                pdir.is_dir()
                and pdir.name not in EXCLUDED_PLUGINS
                and (pdir / "include").is_dir()
            ):
                roots.append(pdir / "include")
    return roots


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
    """Component (and Core sub-module) from a header location path.

    e.g. '.../Core/include/Acts/Surfaces/X.hpp' -> 'Core/Surfaces',
         '.../Plugins/Json/include/ActsPlugins/Json/Y.hpp' -> 'Plugin:Json',
         '.../Fatras/include/ActsFatras/Z.hpp' -> 'Fatras'.
    """
    if not location:
        return "?"
    parts = Path(location).as_posix().split("/")
    if "Plugins" in parts:
        i = parts.index("Plugins")
        if i + 1 < len(parts):
            return "Plugin:" + parts[i + 1]
    if "Fatras" in parts:
        return "Fatras"
    if "Alignment" in parts:
        return "Alignment"
    if "Acts" in parts:  # Core: '.../Core/include/Acts/<module>/file.hpp'
        i = parts.index("Acts")
        return f"Core/{parts[i + 1]}" if i + 1 < len(parts) - 1 else "Core"
    return "?"


def run_doxygen(repo: Path, out: Path, input_dirs: list[Path]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    # Doxygen INPUT accepts a space-separated list of directories.
    doxy_input = " ".join(f'"{d}"' for d in input_dirs)
    env = dict(os.environ, DOXY_OUT=str(out), DOXY_INPUT=doxy_input)
    # Capture output: Doxygen emits many non-fatal cross-reference warnings
    # (our config feeds only headers, so @ref links into the narrative docs do
    # not resolve). They do not affect symbol extraction, so keep them out of
    # the log unless Doxygen actually fails.
    proc = subprocess.run(
        ["doxygen", DOXYFILE], cwd=repo, env=env, capture_output=True, text=True
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise subprocess.CalledProcessError(proc.returncode, "doxygen")
    n = proc.stderr.count("warning:") + proc.stderr.count("error:")
    if n:
        print(f"(doxygen: {n} non-fatal diagnostics suppressed)", file=sys.stderr)


# C++20 permits Unicode identifiers, so a crafted symbol/type name could
# otherwise smuggle bidi-override or zero-width characters into the report
# (Trojan-Source-style visual spoofing) even though every value is rendered
# into a Markdown code span. Strip format (Cf) and control (Cc) code points
# from any text pulled out of Doxygen's XML before it is used as a key or
# rendered.
def _sanitize(text: str) -> str:
    return "".join(c for c in text if unicodedata.category(c) not in ("Cf", "Cc"))


def _name(el, tag: str) -> str:
    return _sanitize(el.findtext(tag) or "")


def _norm_type(el) -> str:
    """Flatten a Doxygen <type>/<param><type> element to normalized text."""
    if el is None:
        return ""
    return _sanitize(" ".join("".join(el.itertext()).split()))


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
    fields: dict[str, str] = {}  # public data member Class::name -> type
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
            name = _name(cd, "compoundname")
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
                # public methods and data members on this documented type
                for md in cd.iter("memberdef"):
                    if md.get("prot") != "public":
                        continue
                    if md.get("kind") == "function":
                        methods += 1
                        callables.update(
                            callable_forms(md, f"{bare}::{_name(md, 'name')}")
                        )
                    elif md.get("kind") == "variable":
                        mname = _name(md, "name")
                        fields[f"{bare}::{mname}"] = _norm_type(md.find("type"))

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
                    mname = _name(md, "name")
                    full = f"{name}::{mname}"
                    if full in ns_member_names:
                        continue
                    ns_member_names.add(full)
                    counts[bucket] += 1
                    symbols.add(f"{bucket} {full}")
                    if mk == "function":
                        callables.update(callable_forms(md, full))
                    mloc = md.find("location")
                    per_module[
                        module_of(mloc.get("file") if mloc is not None else None)
                    ][bucket] += 1

    total = sum(counts.values())
    return {
        "counts": dict(counts),
        "total": total,
        "public_methods": methods,
        "public_fields": len(fields),
        "per_module": {m: dict(c) for m, c in per_module.items()},
        "symbols": sorted(symbols),
        "callables": callables,
        "fields": fields,
    }


BUCKET_ORDER = ["types", "free_functions", "aliases", "variables", "enums", "concepts"]


def render_markdown(data: dict, doxy_version: str | None) -> str:
    c = data["counts"]
    lines = ["## ACTS public API surface", ""]
    lines.append(
        f"**{data['total']}** documented public names in `Acts*::` "
        f"(excluding `detail` / `Experimental`), plus "
        f"**{data['public_methods']}** public methods and "
        f"**{data.get('public_fields', 0)}** public data members on "
        f"documented types."
    )
    scope = data.get("scope")
    if scope:
        lines.append(
            f"\n_Scope: {scope}"
            + (f" via Doxygen {doxy_version}._" if doxy_version else "._")
        )
    elif doxy_version:
        lines.append(f"\n_Doxygen {doxy_version}._")
    lines += ["", "| category | count |", "|---|---:|"]
    for k in BUCKET_ORDER:
        if k in c:
            lines.append(f"| {k.replace('_', ' ')} | {c[k]} |")
    lines.append(f"| **total** | **{data['total']}** |")

    lines += [
        "",
        "<details><summary>By module</summary>",
        "",
        "| module | total | " + " | ".join(BUCKET_ORDER) + " |",
        "|---|---:|" + "|".join(["---:"] * len(BUCKET_ORDER)) + "|",
    ]
    mods = sorted(data["per_module"].items(), key=lambda kv: -sum(kv[1].values()))
    for m, d in mods:
        row = [m, str(sum(d.values()))] + [str(d.get(k, 0)) for k in BUCKET_ORDER]
        lines.append("| " + " | ".join(row) + " |")
    lines += ["", "</details>"]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--repo", default=".", help="repository root (has the Doxyfile)")
    ap.add_argument("--run", action="store_true", help="run doxygen (else use --xml)")
    ap.add_argument(
        "--input", nargs="+", metavar="DIR", help="explicit header root(s) to measure"
    )
    ap.add_argument(
        "--roots-under",
        metavar="DIR",
        help="measure the standard component set (Core, Fatras, Alignment, "
        "and all plugins) resolved under DIR; use for another checkout",
    )
    ap.add_argument("--xml", help="existing Doxygen XML dir to parse")
    ap.add_argument("--json", help="write report JSON here")
    ap.add_argument("--markdown", help="write Markdown here ('-' for stdout)")
    ap.add_argument(
        "--summary",
        action="store_true",
        help="also append the Markdown to $GITHUB_STEP_SUMMARY",
    )
    args = ap.parse_args()

    repo = Path(args.repo).resolve()

    doxy_version = None
    try:
        doxy_version = subprocess.run(
            ["doxygen", "--version"], capture_output=True, text=True
        ).stdout.split()[0]
    except (FileNotFoundError, IndexError):
        pass

    scope = None
    tmp = None
    if args.run:
        if args.input:
            input_dirs = [Path(p).resolve() for p in args.input]
        else:
            input_dirs = standard_roots(
                Path(args.roots_under).resolve() if args.roots_under else repo
            )
        if not input_dirs:
            print("error: no header roots found to measure", file=sys.stderr)
            return 2

        def component(d: Path) -> str:
            parts = d.parts
            if "Plugins" in parts:
                return "Plugin:" + parts[parts.index("Plugins") + 1]
            for comp in ("Core", "Fatras", "Alignment"):
                if comp in parts:
                    return comp
            return d.name

        scope = ", ".join(dict.fromkeys(component(d) for d in input_dirs))
        tmp = Path(tempfile.mkdtemp(prefix="acts-api-surface-"))
        run_doxygen(repo, tmp, input_dirs)
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
    data["scope"] = scope

    md = render_markdown(data, doxy_version)

    if args.json:
        import json

        Path(args.json).write_text(json.dumps(data, indent=2) + "\n")
    if args.markdown == "-" or args.markdown is None:
        print(md)
    elif args.markdown:
        Path(args.markdown).write_text(md)

    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if args.summary and summary:
        with open(summary, "a") as fh:
            fh.write(md)

    return 0


if __name__ == "__main__":
    sys.exit(main())
