#!/usr/bin/env python3
"""Self-tests for check_deprecated_docs.py.

Guards the doc-block association heuristic (the only non-trivial part): walking
up from a ``[[deprecated]]`` past declaration-prefix lines to the Doxygen block,
and stopping at statement boundaries so an unrelated comment above is not
mis-attributed.

Run directly (``python3 CI/public_api/test_check_deprecated_docs.py``); non-zero exit on
failure. Also collectable by pytest.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import check_deprecated_docs as cdd  # noqa: E402


def _reasons(text: str) -> list[str]:
    """Run the checker over a single synthetic header, return violation reasons."""
    with tempfile.TemporaryDirectory() as d:
        repo = Path(d)
        hdr = repo / "Core" / "include" / "Acts" / "T.hpp"
        hdr.parent.mkdir(parents=True)
        hdr.write_text(text)
        vios = cdd.find_violations(repo, ["Core/include"])
        return [v["reason"] for v in vios]


# --- doc_block_for association ---------------------------------------------


def test_paired_simple_block():
    text = (
        "/// Does a thing.\n"
        "/// @deprecated Use bar() instead\n"
        '[[deprecated("Use bar() instead")]] void foo();\n'
    )
    assert _reasons(text) == [], "an adjacent @deprecated block should pair"


def test_template_line_between_doc_and_attribute():
    text = (
        "/// Does a thing.\n"
        "/// @deprecated gone soon\n"
        "template <typename T>\n"
        '[[deprecated("gone soon")]] void foo();\n'
    )
    assert _reasons(text) == [], "template head between block and attr must be skipped"


def test_backslash_command_form_accepted():
    text = "/// \\deprecated old\n[[deprecated]] void foo();\n"
    assert _reasons(text) == []


def test_doc_block_without_command_is_violation():
    text = "/// Does a thing.\n[[deprecated]] void foo();\n"
    assert _reasons(text) == ["no-deprecated-tag"]


def test_no_doc_block_is_violation():
    text = "int x = 0;\n\n[[deprecated]] void foo();\n"
    assert _reasons(text) == ["no-doc-block"]


def test_previous_statement_not_mistaken_for_doc():
    # A @deprecated on the *previous* entity must not satisfy this one.
    text = (
        "/// @deprecated old one\n"
        '[[deprecated("old one")]] void a();\n'
        "\n"
        "[[deprecated]] void b();\n"
    )
    assert _reasons(text) == ["no-doc-block"], "b() is undocumented"


def test_gnu_deprecated_attribute_recognised():
    text = "/// no tag here\n[[gnu::deprecated]] void foo();\n"
    assert _reasons(text) == ["no-deprecated-tag"]


def test_copydoc_block_with_explicit_command_pairs():
    text = (
        "/// @copydoc Base::foo() const\n"
        "/// @deprecated Use the 2D overload instead\n"
        '[[deprecated("Use the 2D overload instead")]] void foo() const;\n'
    )
    assert _reasons(text) == []


def _main() -> int:
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok   {name}")
            except AssertionError as e:
                failures += 1
                print(f"FAIL {name}: {e}")
    print(f"\n{failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())
