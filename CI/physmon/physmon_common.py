import collections
import argparse
import tempfile
import urllib.request
from pathlib import Path

import acts
from acts.examples.odd import getOpenDataDetector

# Provisional Gen3 map, produced locally and hosted as a release asset on a
# personal fork; not yet the ODD repo's official map (which is still
# Gen1-keyed and incompatible with Gen3 surfaces). Hosted externally rather
# than committed: the compressed map is tens of MB, far over this repo's
# check-added-large-files limit.
PROVISIONAL_GEN3_MATERIAL_MAP_URL = (
    "https://github.com/LeaBaumann/acts/releases/download/"
    "gen3-material-map-provisional-v1/odd_gen3_map.json.zst"
)

PhysmonSetup = collections.namedtuple(
    "Setup",
    [
        "detector",
        "trackingGeometry",
        "decorators",
        "field",
        "bFieldInZ",
        "digiConfig",
        "geoSel",
        "gbtsLayerMap",
        "gbtsConnectionTable",
        "outdir",
    ],
)


def makeSetup() -> PhysmonSetup:
    srcdir = Path(__file__).resolve().parent.parent.parent

    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")

    args = parser.parse_args()

    matMapPath = (
        Path(tempfile.gettempdir()) / "odd-material-maps-gen3-provisional.json.zst"
    )
    urllib.request.urlretrieve(PROVISIONAL_GEN3_MATERIAL_MAP_URL, matMapPath)
    matDeco = acts.IMaterialDecorator.fromFile(matMapPath, level=acts.logging.INFO)

    detector = getOpenDataDetector(matDeco)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    field = detector.field
    fieldCache = field.makeCache(acts.MagneticFieldContext())
    bFieldInZ = field.getField(acts.Vector3(0, 0, 0), fieldCache)[2]
    setup = PhysmonSetup(
        detector=detector,
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        digiConfig=srcdir / "Examples/Configs/odd-digi-smearing-config.json",
        geoSel=srcdir / "Examples/Configs/odd-seeding-config.json",
        gbtsLayerMap=srcdir / "Examples/Configs/odd-gbts-layer-map.json",
        gbtsConnectionTable=srcdir / "Examples/Configs/odd-gbts-connection-table.json",
        field=field,
        bFieldInZ=bFieldInZ,
        outdir=Path(args.outdir),
    )

    setup.outdir.mkdir(exist_ok=True)

    return setup
