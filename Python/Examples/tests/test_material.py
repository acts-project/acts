import pytest

import acts

from acts.json import MaterialMapJsonConverter, JsonMaterialDecorator
from acts.examples.odd import getOpenDataDetectorDirectory


@pytest.mark.root
def test_material_root(conf_const):
    with pytest.raises(TypeError):
        acts.root.RootMaterialDecorator()
    fileName = "blubb.root"
    try:
        conf_const(
            acts.root.RootMaterialDecorator,
            level=acts.logging.INFO,
            fileName=fileName,
        )
    except RuntimeError as e:
        assert fileName in str(e)


def test_json_material_decorator():
    config = MaterialMapJsonConverter.Config()
    deco = JsonMaterialDecorator(
        rConfig=config,
        jFileName=str(
            getOpenDataDetectorDirectory() / "config/odd-material-mapping-config.json"
        ),
        level=acts.logging.WARNING,
    )


def test_material_map_export(tmp_path):
    import json
    from acts.examples import GenericDetector
    from acts.json import TrackingGeometryMaterialJsonConverter

    converter = TrackingGeometryMaterialJsonConverter()
    output = tmp_path / "material.json"
    options = TrackingGeometryMaterialJsonConverter.Options()
    options.materialFractionBits = 16
    detector = GenericDetector()
    geometry = detector.trackingGeometry()
    material = acts.TrackingGeometryMaterial.fromGeometry(
        geometry, includeNonMaterial=True
    )
    converter.toFile(material, output, options)
    document = json.loads(output.read_text())
    assert document["header"]["version"] == 1
    assert document["surfaces"]
    assert any(
        entry["material"]["kind"] == "proto-grid" for entry in document["surfaces"]
    )
    material = converter.fromFile(output)
    converter.toFile(material, output, options)
    assert json.loads(output.read_text()) == document

    unsupported = acts.TrackingGeometryMaterial()
    unsupported.volumeMaterials = {acts.GeometryIdentifier(): None}
    with pytest.raises(ValueError, match="surface material only"):
        converter.toFile(unsupported, output)
