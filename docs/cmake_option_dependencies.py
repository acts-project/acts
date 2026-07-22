#!/usr/bin/env python3
"""
Derive the ACTS CMake option auto-enable ("closure") table from a cmake lists
file.

It parses `option(ACTS_... )` declarations and `set_option_if(<target>
<condition>)` implications from the given cmake file, computes for every option
the transitive set of options it auto-enables, and writes a markdown table.

This is the machine-readable counterpart to the `set_option_if` calls in the
top-level `CMakeLists.txt`: enabling an option there can silently switch on
other options during configuration. The table makes that dependency graph
explicit.

Mirrors `docs/parse_cmake_options.py`: use `--write FILE` to update the table in
place between the `<!-- CMAKE_OPTS_DEPS_BEGIN -->` / `<!-- CMAKE_OPTS_DEPS_END
-->` markers, and `--verify` to only check that the file is up to date (used by
the lint job).
"""

import argparse
from pathlib import Path
import re
import sys
import difflib

BEGIN_MARKER = "<!-- CMAKE_OPTS_DEPS_BEGIN -->"
END_MARKER = "<!-- CMAKE_OPTS_DEPS_END -->"


def parse(text, prefix):
    """Return (declared option names, list of (target, condition-tokens))."""
    # drop comments so they cannot leak into parsed conditions
    text = re.sub(r"#.*", "", text)

    options = set(re.findall(rf"option\(\s*({prefix}\w+)", text))

    rules = []
    # set_option_if conditions in this project never contain nested parentheses,
    # so matching up to the first closing paren is sufficient and robust against
    # the one-token-per-line formatting enforced by gersemi.
    for m in re.finditer(r"set_option_if\(([^)]*)\)", text, re.DOTALL):
        tokens = m.group(1).split()
        if not tokens:
            continue
        target, condition = tokens[0], tokens[1:]
        rules.append((target, condition))
    return options, rules


def eval_condition(condition, enabled):
    """Evaluate a set_option_if condition given the set of enabled options.

    The condition is CMake boolean syntax restricted to variable names and the
    AND/OR/NOT operators. CMake and Python share the same operator precedence
    (NOT > AND > OR), so we translate to an equivalent Python expression.
    """
    if not condition:
        return False
    expr = []
    for token in condition:
        if token == "AND":
            expr.append("and")
        elif token == "OR":
            expr.append("or")
        elif token == "NOT":
            expr.append("not")
        elif re.fullmatch(r"[\w:]+", token):
            expr.append("True" if token in enabled else "False")
        else:
            raise ValueError(f"unexpected token in set_option_if condition: {token!r}")
    return bool(eval(" ".join(expr), {"__builtins__": {}}, {}))


def closure(entry, rules):
    """Options transitively auto-enabled when only `entry` is switched on."""
    enabled = {entry}
    changed = True
    while changed:
        changed = False
        for target, condition in rules:
            if target not in enabled and eval_condition(condition, enabled):
                enabled.add(target)
                changed = True
    return sorted(enabled - {entry})


def render_table(rows):
    headers = ("Enabling this option", "…also enables (transitively)")
    widths = [len(h) for h in headers]
    for row in rows:
        for i, col in enumerate(row):
            widths[i] = max(widths[i], len(col))

    def line(cols):
        return (
            "| " + " | ".join(col.ljust(widths[i]) for i, col in enumerate(cols)) + " |"
        )

    out = [line(headers)]
    out.append("|" + "|".join("-" * (w + 2) for w in widths) + "|")
    out += [line(row) for row in rows]
    return "\n".join(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cmakefile", help="Input cmake lists file to parse")
    p.add_argument(
        "--prefix",
        default="ACTS_",
        help="Prefix identifying the options to consider",
    )
    p.add_argument(
        "--write",
        "-w",
        type=Path,
        help="Write table into this file between the CMAKE_OPTS_DEPS markers",
    )
    p.add_argument(
        "--verify",
        "-v",
        action="store_true",
        help="Only verify the target file is up to date, don't write",
    )
    args = p.parse_args()

    options, rules = parse(Path(args.cmakefile).read_text(), args.prefix)

    rows = []
    for option in sorted(options):
        enabled = closure(option, rules)
        if enabled:
            rows.append((option, ", ".join(enabled)))

    output = render_table(rows)

    if not args.write:
        print(output)
        return 0

    if not args.write.exists():
        print(f"Target file {args.write} does not exist")
        return 1

    source = args.write.read_text().split("\n")
    try:
        begin = source.index(BEGIN_MARKER)
        end = source.index(END_MARKER)
    except ValueError:
        print("Markers not found in output file")
        return 1

    if args.verify:
        actual = "\n".join(source[begin + 1 : end])
        if output != actual:
            print("MISMATCH: run docs/cmake_option_dependencies.py --write\n")
            print(
                "\n".join(
                    difflib.unified_diff(
                        actual.split("\n"),
                        output.split("\n"),
                        fromfile="actual",
                        tofile="expected",
                    )
                )
            )
            return 1
        return 0

    out = source[: begin + 1] + output.split("\n") + source[end:]
    args.write.write_text("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
