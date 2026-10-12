"""Path validation and opt-in external nODD integration checks."""

import os
from pathlib import Path
import subprocess
import sys

import pytest

nodd = pytest.importorskip("acts.examples.nodd")


def test_directory_requires_explicit_path(monkeypatch):
    monkeypatch.delenv("NODD_PATH", raising=False)
    with pytest.raises(RuntimeError, match="NODD_PATH"):
        nodd.getNoddDirectory()


def test_directory_override_and_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("NODD_PATH", str(tmp_path / "absent"))
    assert nodd.getNoddDirectory(tmp_path) == tmp_path.resolve()
    with pytest.raises(FileNotFoundError, match="nODD directory"):
        nodd.getNoddDirectory()
    monkeypatch.setenv("NODD_PATH", str(tmp_path))
    assert nodd.getNoddDirectory() == tmp_path.resolve()


def test_missing_compact_fails_before_construction(tmp_path):
    with pytest.raises(FileNotFoundError, match="tracker compact"):
        nodd.getNoddDetector(tmp_path, build_dir="custom-build")


def test_missing_factory_fails_before_xml_parsing(tmp_path):
    compact = tmp_path / "custom.xml"
    compact.write_text("not XML; missing factory must be detected first")
    with pytest.raises(FileNotFoundError, match="factory file"):
        nodd.getNoddDetector(tmp_path, compact_file=compact)


@pytest.fixture
def external_nodd():
    path = os.environ.get("NODD_PATH")
    if not path:
        pytest.skip("Set NODD_PATH to test the external built tracker")
    return Path(path).resolve()


def test_external_detector(external_nodd):
    code = """
from acts.examples.nodd import DD4hepNoddDetector, getNoddDetector
detector = getNoddDetector()
assert isinstance(detector, DD4hepNoddDetector)
assert detector.config.name == 'NoddDetector'
assert len(detector.config.xmlFileNames) == 1
import json
report=json.loads(detector.geometryReport())
assert report['generation']==3
assert len(report['sensors'])>0
assert len({s['geometry_id'] for s in report['sensors']})==len(report['sensors'])
assert {s['ids']['system'] for s in report['sensors']}==set(range(1,8))
print(detector.config.xmlFileNames[0])
"""
    subprocess.run([sys.executable, "-c", code], check=True, timeout=120)


def test_external_material_recording(external_nodd, tmp_path):
    pytest.importorskip("acts.examples.geant4")
    uproot = pytest.importorskip("uproot")
    script = Path(__file__).resolve().parents[1] / "Scripts/nodd_material_recording.py"
    stem = tmp_path / "output" / "material"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "--nodd-dir",
            str(external_nodd),
            "-n",
            "1",
            "-t",
            "8",
            "--seed",
            "228",
            "--material-track-collection",
            "nodd_probe_tracks",
            "-o",
            str(stem),
        ],
        check=True,
        timeout=180,
    )
    with uproot.open(str(stem) + ".root") as output:
        tree = output["nodd_probe_tracks"]
        assert tree.num_entries == 8
        totals = tree["t_X0"].array(library="np")
        assert all(float(value) > 0 and float(value) < float("inf") for value in totals)
