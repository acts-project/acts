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


def test_legacy_material_error():
    path = FIXTURES / "material_map_example.json"
    with pytest.raises(
        ValueError, match="Legacy material format.*ActsMaterialMapMigrate"
    ):
        TrackingGeometryMaterialJsonConverter().fromFile(path)


@pytest.fixture
def material_migrator():
    executable = shutil.which("ActsMaterialMapMigrate")
    if executable is None:
        candidate = (
            Path(acts.__file__).absolute().parents[2] / "bin/ActsMaterialMapMigrate"
        )
        if not candidate.is_file():
            pytest.skip("ActsMaterialMapMigrate is not built")
        executable = str(candidate)

    return executable


def test_material_map_migration(tmp_path, material_migrator):
    executable = material_migrator

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


def test_root_material_map_migration(tmp_path, material_migrator):
    import numpy as np
    import uproot

    source = tmp_path / "material.root"
    output = tmp_path / "material.json"
    geometry_id = (1 << 56) | (256 << 8)
    with uproot.recreate(source) as root:
        root.mktree(
            "HomogeneousMaterial",
            {
                "hGeoId": np.array([geometry_id], dtype=np.int64),
                **{
                    name: np.array([value], dtype=np.float32)
                    for name, value in {
                        "t": 0.25,
                        "x0": 93.7,
                        "l0": 465.2,
                        "A": 28.0855,
                        "Z": 14.0,
                        "rho": 0.00233,
                    }.items()
                },
            },
        )

    def run(path, target):
        return subprocess.run(
            [material_migrator, str(path), str(target)], capture_output=True, text=True
        )

    result = run(source, output)
    if "ROOT input is unavailable" in result.stderr:
        assert result.returncode != 0
        assert "ACTS_BUILD_PLUGIN_ROOT=ON" in result.stderr
        assert not output.exists()
        return
    assert result.returncode == 0, result.stderr
    converter = TrackingGeometryMaterialJsonConverter()
    assert len(converter.fromFile(output).surfaceMaterials) == 1
    document = json.loads(output.read_text())
    entry = document["surfaces"][0]
    assert entry["target"]["geometry_id"] == {"volume": 1, "sensitive": 256}
    assert entry["material"]["slab"]["thickness"] == pytest.approx(0.25)
    assert entry["material"]["slab"]["material"]["radiation_length"] == pytest.approx(
        93.7
    )
    cbor = tmp_path / "material.cbor"
    assert run(source, cbor).returncode == 0
    assert len(converter.fromFile(cbor).surfaceMaterials) == 1

    with uproot.update(source) as root:
        root.mkdir("VolumeMaterial_vol1")
    rejected = run(source, tmp_path / "volume.json")
    assert rejected.returncode != 0
    assert "surface material only" in rejected.stderr
    assert not (tmp_path / "volume.json").exists()

    empty = tmp_path / "empty.root"
    with uproot.recreate(empty):
        pass
    rejected = run(empty, tmp_path / "empty.json")
    assert rejected.returncode != 0
    assert "No material maps" in rejected.stderr
