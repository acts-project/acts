import collections
import argparse
from pathlib import Path

import acts
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

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
    odd_dir = getOpenDataDetectorDirectory()

    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")

    args = parser.parse_args()

    matDeco = acts.IMaterialDecorator.fromFile(
        odd_dir / "data/odd-material-maps.root", level=acts.logging.INFO
    )

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
