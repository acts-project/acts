"""Test cache selection and deletion without accessing GitHub."""

from unittest.mock import Mock, call

import pytest

import prune_ccache


def cache(i, variant="Linux-linux_ubuntu-r2", **changes):
    result = {
        "id": i,
        "key": f"ccache-{variant}-{i:040x}",
        "ref": "refs/heads/main",
        "version": "format-a",
        "size_in_bytes": 100,
        "created_at": f"2026-09-{i:02d}T00:00:00Z",
    }
    return result | changes


def test_only_superseded_main_archives_of_same_variant_and_version():
    old, new = cache(1), cache(2)
    protected = [
        cache(3, ref="refs/pull/123/merge"),
        cache(4, ref="refs/heads/feature"),
        cache(5, version="format-b"),
        cache(6, variant="Linux-linux_ubuntu_extra-r2-clang22-23"),
        cache(7, variant="macOS-macos-r2"),
        cache(8, key="spack-r5-Linux"),
        cache(9, key="ccache-Linux-linux_ubuntu-r2-not-a-sha"),
        cache(10, version=""),
        cache(11, variant="Linux-linux_ubuntu-r3"),
        cache(12, size_in_bytes=0),
    ]
    assert prune_ccache.superseded_caches([new, *protected, old]) == [(new, [old])]
    assert prune_ccache.superseded_caches([old, cache(2, size_in_bytes=0)]) == []


@pytest.mark.parametrize(
    "apply, replacement_present, deleted_ids",
    [(False, True, []), (True, True, [1]), (True, False, [])],
    ids=["dry-run", "replacement-present", "replacement-disappeared"],
)
def test_cache_deletion(monkeypatch, apply, replacement_present, deleted_ids):
    # A concurrent upload must never be included in the deletion plan.
    listing = Mock(return_value=[cache(2), cache(3)] if replacement_present else [])
    delete = Mock()
    monkeypatch.setattr(prune_ccache, "list_caches", listing)
    monkeypatch.setattr(prune_ccache.subprocess, "run", delete)

    prune_ccache.prune("owner/repo", [cache(1), cache(2)], apply=apply)

    assert delete.call_args_list == [
        call(
            ["gh", "api", "--method", "DELETE", f"repos/owner/repo/actions/caches/{i}"],
            check=True,
        )
        for i in deleted_ids
    ]
    if not apply:
        listing.assert_not_called()
