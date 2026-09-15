#!/usr/bin/env python3
"""Exercise the unused-file checker's search semantics using real grep calls."""

import subprocess
import tempfile
import unittest
from pathlib import Path

from check_unused_files import file_can_be_removed


class ReferenceSearchTests(unittest.TestCase):
    def test_references_and_patterns(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "scope with spaces and 'quotes'"
            root.mkdir()
            (root / "references.txt").write_text(
                "#include <used.hpp>\nimport package.module\nfrom helper import f\n"
                "a'b.hpp\n-dash.hpp\n"
            )
            for pattern in (
                "used.hpp",
                r"import .*module",
                "from helper import",
                "a'b.hpp",
                "-dash.hpp",
            ):
                with self.subTest(pattern=pattern):
                    self.assertFalse(file_can_be_removed(pattern, [str(root)]))
            self.assertTrue(file_can_be_removed("unused.hpp", [str(root)]))

    def test_recursive_scopes_and_binary_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first"
            second = root / "second"
            first.mkdir()
            (second / "nested").mkdir(parents=True)
            (first / "binary").write_bytes(b"\0binary_only.hpp\n")
            (second / "nested" / "reference").write_text("nested.hpp\n")
            scopes = [str(first), str(second)]
            self.assertFalse(file_can_be_removed("nested.hpp", scopes))
            self.assertTrue(file_can_be_removed("binary_only.hpp", scopes))
            self.assertTrue(file_can_be_removed("nested.hpp", [str(first)]))

    def test_search_errors_fail_the_check(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(subprocess.CalledProcessError):
                file_can_be_removed("unused.hpp", [str(Path(directory) / "missing")])


if __name__ == "__main__":
    unittest.main()
