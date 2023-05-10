import collections

import acts

Setup = collections.namedtuple("Setup", ["detector", "trackingGeometry", "decorators", "field", "digiConfig", "geoSel", "outdir"])


def makeSetup() -> Setup:
    u = acts.UnitConstants

    setup = Setup()
    setup.matDeco = acts.IMaterialDecorator.fromFile(
        srcdir / "thirdparty/OpenDataDetector/data/odd-material-maps.root",
        level=acts.logging.INFO,
    )
    setup.detector, setup.trackingGeometry, setup.decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory(), setup.matDeco
    )
    setup.digiConfig = srcdir / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json"
    setup.geoSel = srcdir / "thirdparty/OpenDataDetector/config/odd-seeding-config.json"

    setup.field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")

    args = parser.parse_args()

    setup.outdir = Path(args.outdir)
    setup.outdir.mkdir(exist_ok=True)

    return setup


