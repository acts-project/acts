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
