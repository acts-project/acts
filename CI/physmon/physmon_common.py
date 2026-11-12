import collections
import argparse
import shutil
from pathlib import Path

import acts
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

u = acts.UnitConstants

SRC_DIR = Path(__file__).resolve().parent.parent.parent
GBTS_LAYER_MAP = SRC_DIR / "Examples/Configs/odd-gbts-layer-map.json"
GBTS_CONNECTION_TABLE = SRC_DIR / "Examples/Configs/odd-gbts-connection-table.json"

# Starting uncertainties of the estimated seed parameters, the same in every
# trackfinding chain.
SEED_PARAMETER_SIGMAS = dict(
    initialSigmas=[
        1 * u.mm,
        1 * u.mm,
        1 * u.degree,
        1 * u.degree,
        0 * u.e / u.GeV,
        1 * u.ns,
    ],
    initialSigmaQoverPt=0.1 * u.e / u.GeV,
    initialSigmaPtRel=0.1,
    initialVarInflation=[1.0] * 6,
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
        "outdir",
        "threads",
        "variant",
    ],
)


def makeSetup(variants=None) -> PhysmonSetup:
    """Parse the command line and set up the detector.

    `variants` are the `SeedingAlgorithm` member names a trackfinding script
    supports; they enable `--variant`, see seedingVariants().
    """
    odd_dir = getOpenDataDetectorDirectory()

    parser = argparse.ArgumentParser()
    parser.add_argument("outdir")
    parser.add_argument(
        "--threads",
        type=int,
        default=-1,
        help="Number of threads Snakemake can use, -1 for all",
    )
    if variants is not None:
        parser.add_argument(
            "--variant",
            choices=variants,
            default=None,
            help="Run only this seeding variant and write its files directly to outdir.",
        )

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
        digiConfig=SRC_DIR / "Examples/Configs/odd-digi-smearing-config.json",
        geoSel=SRC_DIR / "Examples/Configs/odd-seeding-config.json",
        field=field,
        bFieldInZ=bFieldInZ,
        outdir=Path(args.outdir),
        threads=args.threads,
        variant=getattr(args, "variant", None),
    )

    setup.outdir.mkdir(exist_ok=True)

    return setup


def seedingVariants(setup: PhysmonSetup, variants):
    """The (SeedingAlgorithm, output directory) pairs a trackfinding script runs.

    With `--variant` only that one runs, with its files directly in outdir:
    snakemake runs one job per variant, each with its own outdir. Without it
    every variant runs, one after the other, into outdir/<variant>/, which is
    what a manual invocation outside snakemake wants.
    """
    from acts.examples.reconstruction import SeedingAlgorithm

    if setup.variant is not None:
        return [(getattr(SeedingAlgorithm, setup.variant), setup.outdir)]
    return [(getattr(SeedingAlgorithm, v), setup.outdir / v) for v in variants]


def copyOutputs(outputs, source: Path, outdir: Path):
    """Copy the files a chain wrote to `source` into outdir.

    `outputs` maps each file, relative to `source`, to the name it gets.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    for written, name in outputs.items():
        assert (source / written).exists(), f"Output not found: {source / written}"
        shutil.copy(source / written, outdir / name)
