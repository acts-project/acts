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
        "digiConfig",
        "geoSel",
        "gbtsLayerMap",
        "gbtsConnectionTable",
        "outdir",
    ],
)


def makeSetup() -> PhysmonSetup:
    u = acts.UnitConstants
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
    setup = PhysmonSetup(
        detector=detector,
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        digiConfig=srcdir / "Examples/Configs/odd-digi-smearing-config.json",
        geoSel=srcdir / "Examples/Configs/odd-seeding-config.json",
        gbtsLayerMap=srcdir / "Examples/Configs/odd-gbts-layer-map.json",
        gbtsConnectionTable=srcdir / "Examples/Configs/odd-gbts-connection-table.json",
        field=acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T)),
        outdir=Path(args.outdir),
    )

    setup.outdir.mkdir(exist_ok=True)

    return setup
