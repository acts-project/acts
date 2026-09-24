"""Versioned material files and the explicit legacy compatibility bindings."""

import json
import shutil
import subprocess
from pathlib import Path

import pytest
import acts

acts_json = pytest.importorskip("acts.json")
TrackingGeometryMaterialJsonConverter = acts_json.TrackingGeometryMaterialJsonConverter
MaterialMapDecorator = acts_json.MaterialMapDecorator

FIXTURES = Path(__file__).resolve().parents[3] / "docs/examples"


def test_material_document(tmp_path):
    converter = TrackingGeometryMaterialJsonConverter()
    material = converter.fromFile(FIXTURES / "material-map-v1/minimal.json")
    material.description = "Python round trip"
    options = TrackingGeometryMaterialJsonConverter.Options()
    options.materialFractionBits = 16
    for suffix in (".json", ".cbor"):
        path = tmp_path / ("material" + suffix)
        converter.toFile(material, path, options)
        assert converter.fromFile(path).description == "Python round trip"
        assert isinstance(acts.IMaterialDecorator.fromFile(path), MaterialMapDecorator)
    encoded = json.loads((tmp_path / "material.json").read_text())
    assert encoded["header"]["version"] == 1
    assert encoded["header"]["description"] == "Python round trip"


def test_legacy_material_warning_and_error():
    from acts.json import JsonMaterialDecorator, MaterialMapJsonConverter

    path = FIXTURES / "material_map_example.json"
    with pytest.raises(
        ValueError, match="Legacy material format.*ActsMaterialMapMigrate"
    ):
        TrackingGeometryMaterialJsonConverter().fromFile(path)
    config = MaterialMapJsonConverter.Config()
    with pytest.warns(DeprecationWarning, match="Legacy JSON material APIs"):
        MaterialMapJsonConverter(config, acts.logging.WARNING)
    with pytest.warns(DeprecationWarning, match="Legacy JSON material APIs"):
        JsonMaterialDecorator(config, str(path), acts.logging.WARNING)
    with pytest.warns(DeprecationWarning, match="Legacy JSON material APIs"):
        acts.IMaterialDecorator.fromFile(path)


def test_material_map_migration(tmp_path):
    executable = shutil.which("ActsMaterialMapMigrate")
    if executable is None:
        candidate = (
            Path(acts.__file__).absolute().parents[2] / "bin/ActsMaterialMapMigrate"
        )
        if not candidate.is_file():
            pytest.skip("ActsMaterialMapMigrate is not built")
        executable = str(candidate)

    def run(*args):
        return subprocess.run(
            [executable, *map(str, args)], capture_output=True, text=True
        )

    original = FIXTURES / "material_map_example.json"
    output = tmp_path / "new.json"
    rejected = run(original, output)
    assert rejected.returncode != 0
    assert "surface material only" in rejected.stderr
    assert not output.exists()
    legacy = json.loads(original.read_text())
    legacy["Volumes"]["entries"] = []
    legacy["KeyedSurfaces"] = [
        {
            "key": "migration-test",
            "geometry_id": 0,
            "material": legacy["Surfaces"]["entries"][0]["value"]["material"],
        }
    ]
    source = tmp_path / "legacy.json"
    source.write_text(json.dumps(legacy))
    migrated = run(source, output)
    assert migrated.returncode == 0, migrated.stderr
    material = TrackingGeometryMaterialJsonConverter().fromFile(output)
    assert "migration-test" in material.keyedSurfaces
    assert len(material.surfaceMaterials) == len(legacy["Surfaces"]["entries"])
    quantized = tmp_path / "quantized.json"
    result = run(source, quantized, "--material-fraction-bits", 16)
    assert result.returncode == 0, result.stderr
    TrackingGeometryMaterialJsonConverter().fromFile(quantized)
    for args in [
        (output, tmp_path / "invalid.json"),
        (source, source),
        (source, output, "--material-fraction-bits", "24"),
    ]:
        assert run(*args).returncode != 0
