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
        "threads",
    ],
)


def makeSetup() -> PhysmonSetup:
    srcdir = Path(__file__).resolve().parent.parent.parent
    odd_dir = getOpenDataDetectorDirectory()

    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    # -1 keeps the legacy "use all hardware threads" behavior for manual runs;
    # the Snakefile passes a concrete value matching the rule's thread budget.
    parser.add_argument("--threads", type=int, default=-1)

    # parse_known_args so callers can layer their own argparse on top
    # (e.g. physmon_trackfinding_1muon.py adds a seeding-variant argument)
    args, _ = parser.parse_known_args()

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
        threads=args.threads,
    )

    setup.outdir.mkdir(exist_ok=True)

    return setup
