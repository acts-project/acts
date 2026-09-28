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
    monkeypatch.setattr(snapshot.subprocess, "run", lambda *args, **kwargs: None)
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
    "event,ref", [("pull_request", "refs/pull/1/merge"), ("push", "refs/heads/feature")]
)
def test_only_main_pushes_publish(monkeypatch, event, ref):
    monkeypatch.setenv("GITHUB_EVENT_NAME", event)
    monkeypatch.setenv("GITHUB_REF", ref)
    monkeypatch.setattr(sys, "argv", ["snapshot", "publish"])
    monkeypatch.setattr(
        snapshot, "make_client", lambda *_: pytest.fail("Must not connect")
    )
    assert snapshot.main() == 0


def test_cache_failure_does_not_fail_build(monkeypatch):
    def unavailable(*args):
        raise OSError("Unavailable")

    monkeypatch.setattr(snapshot, "make_client", unavailable)
    monkeypatch.setattr(sys, "argv", ["snapshot", "restore"])
    assert snapshot.main() == 0
