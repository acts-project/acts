"""Exercise snapshot replacement and restoration without remote credentials."""

import io
import json
from pathlib import Path
import sys
import tarfile

import pytest

import ccache_snapshot as snapshot


class MissingObject(Exception):
    response = {"Error": {"Code": "NoSuchKey"}}


class FakeS3:
    def __init__(self):
        self.objects = {}
        self.writes = []
        self.fail_upload = False

    def get_paginator(self, operation):
        assert operation == "list_objects_v2"
        return self

    def paginate(self, *, Bucket, Prefix):
        objects = [
            {"Key": key, "Size": len(value)}
            for key, value in self.objects.items()
            if key.startswith(Prefix)
        ]
        # Deliberately small pages exercise pagination with ordinary fixtures.
        for offset in range(0, len(objects), 2):
            yield {"Contents": objects[offset : offset + 2]}

    def delete_objects(self, *, Bucket, Delete):
        for item in Delete["Objects"]:
            del self.objects[item["Key"]]
        return {"Deleted": Delete["Objects"]}

    def get_object(self, *, Bucket, Key):
        if Key not in self.objects:
            raise MissingObject()
        return {"Body": io.BytesIO(self.objects[Key])}

    def download_file(self, bucket, key, filename):
        Path(filename).write_bytes(self.objects[key])

    def upload_file(self, filename, bucket, key):
        if self.fail_upload:
            raise OSError("Simulated upload failure")
        self.writes.append(key)
        self.objects[key] = Path(filename).read_bytes()

    def put_object(self, *, Bucket, Key, Body, **kwargs):
        self.writes.append(Key)
        self.objects[Key] = Body


@pytest.fixture
def client(monkeypatch):
    # Transfer must work on hosts without ccache or other build tools.
    monkeypatch.setenv("PATH", "")
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)
    return FakeS3()


@pytest.fixture
def cache(tmp_path):
    directory = tmp_path / "cache"
    (directory / "a" / "b").mkdir(parents=True)
    (directory / "a" / "b" / "manifestM").write_bytes(b"manifest")
    (directory / "a" / "b" / "resultR").write_bytes(b"compiled output")
    return directory


def publish(client, cache, run=10, attempt=1):
    snapshot.publish(client, "cache", "pilot/", cache, run, 1234, attempt)


def test_round_trip_preserves_manifests_and_results(client, cache, tmp_path):
    publish(client, cache)
    manifest = snapshot.read_manifest(client, "cache", "pilot/")
    assert client.writes == [manifest["key"], "pilot/latest.json"]
    restored = tmp_path / "restored"
    snapshot.restore(client, "cache", "pilot/", restored)
    for filename in ("manifestM", "resultR"):
        assert (restored / "a/b" / filename).read_bytes() == (
            cache / "a/b" / filename
        ).read_bytes()


def test_missing_snapshot_starts_empty(client, tmp_path):
    restored = tmp_path / "restored"
    snapshot.restore(client, "cache", "pilot/", restored)
    assert list(restored.iterdir()) == []


def test_checksum_failure_leaves_cache_empty(client, cache, tmp_path):
    publish(client, cache)
    manifest = snapshot.read_manifest(client, "cache", "pilot/")
    client.objects[manifest["key"]] += b"corruption"
    restored = tmp_path / "restored"
    with pytest.raises(ValueError, match="checksum"):
        snapshot.restore(client, "cache", "pilot/", restored)
    assert list(restored.iterdir()) == []


@pytest.mark.parametrize("entry", ["../escape", "/absolute", "symlink", "hardlink"])
def test_unsafe_archive_is_never_installed(client, cache, tmp_path, entry):
    publish(client, cache)
    archive = tmp_path / "bad.tar"
    with tarfile.open(archive, "w:") as output:
        output.addfile(tarfile.TarInfo("safe-first"), io.BytesIO())
        member = tarfile.TarInfo(entry)
        if entry in ("symlink", "hardlink"):
            member.type = tarfile.SYMTYPE if entry == "symlink" else tarfile.LNKTYPE
            member.linkname = "../escape"
        output.addfile(member, io.BytesIO())
    manifest = snapshot.read_manifest(client, "cache", "pilot/")
    client.objects[manifest["key"]] = archive.read_bytes()
    manifest["sha256"] = snapshot.digest(archive)
    client.objects["pilot/latest.json"] = json.dumps(manifest).encode()
    restored = tmp_path / "restored"
    with pytest.raises(ValueError, match="Unsafe"):
        snapshot.restore(client, "cache", "pilot/", restored)
    assert list(restored.iterdir()) == []
    assert not (tmp_path / "escape").exists()


def test_restore_does_not_discard_existing_cache(client, cache):
    with pytest.raises(ValueError, match="empty"):
        snapshot.restore(client, "cache", "pilot/", cache)
    assert (cache / "a/b/manifestM").read_bytes() == b"manifest"


def test_failed_upload_preserves_pointer(client, cache):
    publish(client, cache)
    previous = client.objects["pilot/latest.json"]
    client.fail_upload = True
    with pytest.raises(OSError):
        publish(client, cache, run=11)
    assert client.objects["pilot/latest.json"] == previous


@pytest.mark.parametrize("run,attempt", [(9, 5), (10, 1)])
def test_old_or_duplicate_run_skips_already_published_snapshot(
    client, cache, run, attempt
):
    publish(client, cache)
    writes = client.writes.copy()
    publish(client, cache, run, attempt)
    assert client.writes == writes


def test_new_attempt_replaces_snapshot(client, cache):
    publish(client, cache)
    previous = snapshot.read_manifest(client, "cache", "pilot/")
    publish(client, cache, attempt=2)
    current = snapshot.read_manifest(client, "cache", "pilot/")
    assert current["run_attempt"] == 2
    assert current["key"] != previous["key"]
    assert previous["key"] in client.objects


def test_access_denied_is_not_treated_as_missing(client, monkeypatch):
    class Denied(Exception):
        response = {"Error": {"Code": "AccessDenied"}}

    def denied(**kwargs):
        raise Denied()

    monkeypatch.setattr(client, "get_object", denied)
    with pytest.raises(Denied):
        snapshot.read_manifest(client, "cache", "pilot/")


def test_manifest_cannot_reference_another_namespace(client, cache):
    publish(client, cache)
    manifest = snapshot.read_manifest(client, "cache", "pilot/")
    manifest["key"] = manifest["key"].replace("pilot/", "another/")
    client.objects["pilot/latest.json"] = json.dumps(manifest).encode()
    with pytest.raises(ValueError, match="Invalid"):
        snapshot.read_manifest(client, "cache", "pilot/")


@pytest.mark.parametrize(
    "event,ref,opt_in,allowed",
    [
        ("push", "refs/heads/main", "false", True),
        ("push", "refs/heads/feature", "true", False),
        ("pull_request", "refs/heads/main", "true", False),
        ("merge_group", "refs/heads/main", "true", False),
        ("schedule", "refs/heads/main", "false", False),
        ("schedule", "refs/heads/main", "true", True),
        ("workflow_dispatch", "refs/heads/main", "false", False),
        ("workflow_dispatch", "refs/heads/main", "true", True),
        ("workflow_dispatch", "refs/heads/feature", "true", False),
    ],
)
def test_publication_policy(monkeypatch, event, ref, opt_in, allowed):
    monkeypatch.setenv("GITHUB_EVENT_NAME", event)
    monkeypatch.setenv("GITHUB_REF", ref)
    monkeypatch.setenv("CACHE_ALLOW_MAIN_NON_PUSH", opt_in)
    monkeypatch.setattr(sys, "argv", ["snapshot", "publish"])
    for key, value in {
        "CACHE_BUCKET": "cache",
        "CACHE_PREFIX": "pilot/",
        "CCACHE_DIR": "/unused",
        "GITHUB_RUN_NUMBER": "1",
        "GITHUB_RUN_ID": "1",
        "GITHUB_RUN_ATTEMPT": "1",
    }.items():
        monkeypatch.setenv(key, value)
    connected, published = [], []
    monkeypatch.setattr(
        snapshot, "make_client", lambda writable: connected.append(writable)
    )
    monkeypatch.setattr(snapshot, "publish", lambda *args: published.append(args))
    assert snapshot.main() == 0
    assert bool(connected) == allowed
    assert bool(published) == allowed


def test_cache_failure_does_not_fail_build(monkeypatch):
    def unavailable(*args):
        raise OSError("Unavailable")

    monkeypatch.setattr(snapshot, "make_client", unavailable)
    monkeypatch.setattr(sys, "argv", ["snapshot", "restore"])
    assert snapshot.main() == 0


ROOT = "acts-sccache/ccache-snapshots/acts-project/acts/"
VARIANT = ROOT + "builds/linux/v1/"


def seed_snapshots(client, cache, variant=VARIANT):
    snapshot.publish(client, "cache", variant, cache, 1, 1, 1)
    old = snapshot.read_manifest(client, "cache", variant)["key"]
    snapshot.publish(client, "cache", variant, cache, 2, 2, 1)
    latest = snapshot.read_manifest(client, "cache", variant)["key"]
    return old, latest


@pytest.mark.parametrize("dry_run", [False, True])
def test_prune_preserves_latest_and_unrelated_objects(client, cache, dry_run):
    old, latest = seed_snapshots(client, cache)
    other_old, other_latest = seed_snapshots(client, cache, ROOT + "macos/v1/")
    client.objects[ROOT + "builds/linux/v1/snapshots/unrecognized.tar"] = b"keep"
    client.objects["another-repository/snapshots/archive.tar"] = b"keep"
    before = client.objects.copy()
    snapshot.prune(client, "cache", ROOT, dry_run)
    expected = before.copy()
    if not dry_run:
        del expected[old]
        del expected[other_old]
    assert client.objects == expected
    assert latest in client.objects and other_latest in client.objects


@pytest.mark.parametrize("manifest", [None, b"invalid json"])
def test_prune_missing_or_invalid_pointer_preserves_variant(client, cache, manifest):
    seed_snapshots(client, cache)
    if manifest is None:
        del client.objects[VARIANT + "latest.json"]
    else:
        client.objects[VARIANT + "latest.json"] = manifest
    before = client.objects.copy()
    with pytest.raises(RuntimeError):
        snapshot.prune(client, "cache", ROOT)
    assert client.objects == before


def test_prune_rechecks_pointer_and_excludes_later_uploads(client, cache, monkeypatch):
    old, latest = seed_snapshots(client, cache)
    previous = json.loads(client.objects[VARIANT + "latest.json"])
    previous["key"] = old
    get = client.get_object
    reads = 0

    def concurrent_publish(**kwargs):
        nonlocal reads
        reads += 1
        if reads == 1:
            # A completed upload after listing must not become a candidate.
            client.objects[VARIANT + "snapshots/3-1-" + "a" * 32 + ".tar"] = b"new"
        if reads == 2:
            # Overlapping publishers can move the pointer backwards, too.
            client.objects[VARIANT + "latest.json"] = json.dumps(previous).encode()
        return get(**kwargs)

    monkeypatch.setattr(client, "get_object", concurrent_publish)
    before = set(client.objects)
    snapshot.prune(client, "cache", ROOT)
    assert before <= set(client.objects)
    assert len(client.objects) == len(before) + 1
    assert old in client.objects and latest in client.objects


def test_prune_batches_deletions(client, cache, monkeypatch):
    _, latest = seed_snapshots(client, cache)
    for number in range(1005):
        client.objects[VARIANT + f"snapshots/3-1-{number:032x}.tar"] = b"old"
    delete = client.delete_objects
    batches = []

    def limited_delete(**kwargs):
        batches.append(len(kwargs["Delete"]["Objects"]))
        assert batches[-1] <= 1000
        return delete(**kwargs)

    monkeypatch.setattr(client, "delete_objects", limited_delete)
    snapshot.prune(client, "cache", ROOT)
    assert len(batches) == 2
    assert set(client.objects) == {latest, VARIANT + "latest.json"}


def test_prune_partial_deletion_is_failure(client, cache, monkeypatch, capsys):
    old, latest = seed_snapshots(client, cache)

    def denied(**kwargs):
        return {"Errors": [{"Key": old, "Code": "AccessDenied"}]}

    monkeypatch.setattr(client, "delete_objects", denied)
    with pytest.raises(RuntimeError):
        snapshot.prune(client, "cache", ROOT)
    assert old in client.objects and latest in client.objects
    assert "deleted 0 archives" in capsys.readouterr().out


@pytest.mark.parametrize(
    "repository,ref,event,prefix,allowed",
    [
        ("acts-project/acts", "refs/heads/main", "schedule", ROOT, True),
        ("acts-project/acts", "refs/heads/main", "workflow_dispatch", ROOT, True),
        ("acts-project/acts", "refs/heads/main", "push", ROOT, False),
        ("acts-project/acts", "refs/heads/main", "pull_request", ROOT, False),
        ("acts-project/acts", "refs/heads/feature", "workflow_dispatch", ROOT, False),
        ("fork/acts", "refs/heads/main", "schedule", ROOT, False),
        ("acts-project/acts", "refs/heads/main", "schedule", "acts-sccache/", False),
    ],
)
def test_prune_policy(monkeypatch, repository, ref, event, prefix, allowed):
    for key, value in {
        "GITHUB_REPOSITORY": repository,
        "GITHUB_REF": ref,
        "GITHUB_EVENT_NAME": event,
        "CACHE_PREFIX": prefix,
        "CACHE_BUCKET": "cache",
    }.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(sys, "argv", ["snapshot", "prune", "--dry-run"])
    connected, pruned = [], []
    monkeypatch.setattr(
        snapshot, "make_client", lambda writable: connected.append(writable)
    )
    monkeypatch.setattr(snapshot, "prune", lambda *args: pruned.append(args))
    assert snapshot.main() == (0 if allowed else 1)
    assert bool(connected) == allowed
    assert bool(pruned) == allowed
    if allowed:
        assert pruned[0][-1] is True


def test_prune_failure_fails_job(monkeypatch):
    for key, value in {
        "GITHUB_REPOSITORY": "acts-project/acts",
        "GITHUB_REF": "refs/heads/main",
        "GITHUB_EVENT_NAME": "schedule",
        "CACHE_PREFIX": ROOT,
    }.items():
        monkeypatch.setenv(key, value)

    def unavailable(*args):
        raise OSError("Unavailable")

    monkeypatch.setattr(snapshot, "make_client", unavailable)
    monkeypatch.setattr(sys, "argv", ["snapshot", "prune"])
    assert snapshot.main() == 1
