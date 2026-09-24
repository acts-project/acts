"""Versioned material files and the explicit legacy compatibility bindings."""

import json
from pathlib import Path

import pytest
import acts

acts_json = pytest.importorskip("acts.json")
TrackingGeometryMaterialJsonConverter = acts_json.TrackingGeometryMaterialJsonConverter

FIXTURES = Path(__file__).resolve().parents[3] / "docs/examples"


def test_material_document(tmp_path):
    converter = TrackingGeometryMaterialJsonConverter()
    material = converter.fromFile(FIXTURES / "material-map-v1/minimal.json")
    material.description = "Python round trip"
    options = TrackingGeometryMaterialJsonConverter.Options()
    for suffix in (".json", ".cbor"):
        path = tmp_path / ("material" + suffix)
        converter.toFile(material, path, options)
        assert converter.fromFile(path).description == "Python round trip"
    encoded = json.loads((tmp_path / "material.json").read_text())
    assert encoded["header"]["version"] == 1
    assert encoded["header"]["description"] == "Python round trip"
