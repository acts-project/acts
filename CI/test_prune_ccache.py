#!/usr/bin/env python3
"""Test retention boundaries and deletion behavior without accessing GitHub."""

from contextlib import redirect_stdout
import io
import json
import subprocess
import unittest
from unittest.mock import patch

from prune_ccache import list_caches, prune, superseded_caches


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


class RetentionTests(unittest.TestCase):
    def test_only_superseded_main_archives_of_same_variant_and_version(self):
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
        ]
        self.assertEqual(superseded_caches([new, *protected, old]), [(new, [old])])

    def test_empty_replacement_does_not_displace_last_cache(self):
        self.assertEqual(superseded_caches([cache(1), cache(2, size_in_bytes=0)]), [])

    def test_order_by_creation_not_input_order_or_last_access(self):
        old, new = cache(1, last_accessed_at="2026-09-20"), cache(2)
        self.assertEqual(superseded_caches([new, old]), [(new, [old])])

    def test_dry_run_does_not_call_github(self):
        with patch("prune_ccache.list_caches") as listing, patch(
            "prune_ccache.subprocess.run"
        ) as delete, redirect_stdout(io.StringIO()):
            self.assertEqual(prune("owner/repo", [cache(1), cache(2)]), 100)
        listing.assert_not_called()
        delete.assert_not_called()

    def test_disappearing_replacement_preserves_old_cache(self):
        with patch("prune_ccache.list_caches", return_value=[]), patch(
            "prune_ccache.subprocess.run"
        ) as delete, redirect_stdout(io.StringIO()):
            self.assertEqual(prune("owner/repo", [cache(1), cache(2)], True), 0)
        delete.assert_not_called()

    def test_delete_only_old_ids_from_snapshot(self):
        with patch(
            "prune_ccache.list_caches", return_value=[cache(2), cache(3)]
        ), patch("prune_ccache.subprocess.run") as delete, redirect_stdout(
            io.StringIO()
        ):
            self.assertEqual(prune("owner/repo", [cache(1), cache(2)], True), 100)
        self.assertEqual(
            delete.call_args.args[0][-1], "repos/owner/repo/actions/caches/1"
        )
        self.assertEqual(delete.call_count, 1)

    def test_all_pages_are_combined_before_planning(self):
        pages = [{"actions_caches": [cache(1)]}, {"actions_caches": [cache(2)]}]
        with patch(
            "prune_ccache.subprocess.check_output", return_value=json.dumps(pages)
        ):
            self.assertEqual(list_caches("owner/repo"), [cache(1), cache(2)])

    def test_api_error_does_not_become_an_empty_listing(self):
        with patch(
            "prune_ccache.subprocess.check_output",
            side_effect=subprocess.CalledProcessError(1, "gh"),
        ):
            with self.assertRaises(subprocess.CalledProcessError):
                list_caches("owner/repo")


if __name__ == "__main__":
    unittest.main()
