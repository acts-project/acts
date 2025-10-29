import pytest
from pathlib import Path

import acts

from helpers import jsonEnabled


@pytest.mark.root
def test_material_root(conf_const):
    from acts import root

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


@pytest.mark.skipif(not jsonEnabled, reason="JSON not set up")
def test_json_material_decorator():
    import acts.json

    config = acts.json.MaterialMapJsonConverter.Config()
    deco = acts.json.JsonMaterialDecorator(
        rConfig=config,
        jFileName=str(
            Path(__file__).parent.parent.parent.parent
            / "thirdparty/OpenDataDetector/config/odd-material-mapping-config.json"
        ),
        level=acts.logging.WARNING,
    )
