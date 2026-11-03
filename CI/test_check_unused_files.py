"""Exercise reference matching and search failures with real grep calls."""

import subprocess

import pytest

from check_unused_files import file_can_be_removed


def test_reference_matching(tmp_path):
    (tmp_path / "references.txt").write_text(
        "#include <used.hpp>\nimport package.module\n"
    )
    scope = [str(tmp_path)]
    assert not file_can_be_removed("used.hpp", scope)
    assert not file_can_be_removed(r"import .*module", scope)
    assert file_can_be_removed("unused.hpp", scope)


def test_search_failure(tmp_path):
    with pytest.raises(subprocess.CalledProcessError):
        file_can_be_removed("unused.hpp", [str(tmp_path / "missing")])
