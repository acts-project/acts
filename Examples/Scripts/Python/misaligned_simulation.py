#!/usr/bin/env python3
"""Simulate a detector that is placed differently from what reconstruction assumes.

The `AlignmentDecorator` writes its alignment into the simulation geometry context
only, so Fatras and the digitization see a shifted layer while seeding, the CKF and
the fit stay on the nominal geometry. The shift then shows up as a bias in the
track state residuals of that layer, which is what an alignment procedure has to
recover.

Writing measurements out of a simulation job and reading them back into a separate
reconstruction job is not necessary for this - both geometries live in the same
sequence.
"""

import argparse
import os
from pathlib import Path

import acts
import acts.examples
from acts import UnitConstants as u
from acts.examples import (
    GaussianVertexGenerator,
    RandomNumbers,
    Sequencer,
    StructureSelector,
    TelescopeDetector,
)
from acts.examples.alignment import (
    AlignmentDecorator,
    AlignmentGeneratorGlobalShift,
    GeoIdAlignmentStore,
)
from acts.examples.reconstruction import (
    CkfConfig,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    TrackSelectorConfig,
    addCKFTracks,
    addSeeding,
)
from acts.examples.simulation import (
    EtaConfig,
    MomentumConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    PhiConfig,
    addDigiParticleSelection,
    addDigitization,
    addFatras,
    addParticleGun,
)

# The telescope layout of millepede_alignment.py: square sensors in the x-z plane,
# geo id layer running 2, 4, ..., 18 for the nine layers, all in volume 1.
LAYER_POSITIONS = [30, 60, 90, 120, 150, 180, 210, 240, 270]


def addMisalignment(
    s: Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    layer: int,
    shift: acts.Vector3,
    target: AlignmentDecorator.Target,
    logLevel: acts.logging.Level = acts.logging.WARNING,
) -> AlignmentDecorator:
    """Shift one layer, and decorate the selected geometry context(s) with it."""
    geoId = acts.GeometryIdentifier(volume=1, layer=layer, sensitive=1)

    generator = AlignmentGeneratorGlobalShift()
    generator.shift = shift

    cfg = AlignmentDecorator.Config()
    cfg.target = target
    cfg.nominalStore = GeoIdAlignmentStore(
        StructureSelector(trackingGeometry).selectedTransforms(
            acts.GeometryContext.dangerouslyDefaultConstruct(), geoId
        )
    )
    # A single IOV covering the whole run - the point here is the sim/reco split,
    # not time dependence.
    cfg.iovGenerators = [((0, 10000000), generator)]

    decorator = AlignmentDecorator(cfg, logLevel)
    s.addContextDecorator(decorator)
    return decorator


def addTelescopeChain(
    s: Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    outputDir: Path,
    rnd: RandomNumbers,
) -> None:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    addParticleGun(
        s,
        MomentumConfig(10 * u.GeV, 100 * u.GeV, transverse=True),
        EtaConfig(-0.3, 0.3),
        # the telescope points along +y
        PhiConfig(60 * u.degree, 120 * u.degree),
        ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
        vtxGen=GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(5.0 * u.mm, 0.0 * u.mm, 5.0 * u.mm, 0.0 * u.ns),
        ),
        multiplicity=1,
        rnd=rnd,
    )
    addFatras(
        s,
        trackingGeometry,
        field,
        enableInteractions=True,
        rnd=rnd,
    )
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=srcdir / "Examples/Configs/telescope-digi-smearing-config.json",
        rnd=rnd,
    )
    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(measurements=(3, None), removeNeutral=True),
    )
    addSeeding(
        s,
        trackingGeometry,
        field,
        seedFinderConfigArg=SeedFinderConfigArg(
            r=(20 * u.mm, 200 * u.mm),
            deltaR=(1 * u.mm, 300 * u.mm),
            collisionRegion=(-250 * u.mm, 250 * u.mm),
            z=(-100 * u.mm, 100 * u.mm),
            maxSeedsPerSpM=1,
            sigmaScattering=5,
            radLengthPerSeed=0.1,
            minPt=0.5 * u.GeV,
            impactMax=3 * u.mm,
        ),
        seedFinderOptionsArg=SeedFinderOptionsArg(bFieldInZ=2 * u.T),
        seedingAlgorithm=SeedingAlgorithm.GridTriplet,
        initialSigmas=[
            3 * u.mm,
            3 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0 * u.e / u.GeV,
            1 * u.ns,
        ],
        initialSigmaQoverPt=0.1 * u.e / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1.0] * 6,
        geoSelectionConfigFile=srcdir
        / "Examples/Configs/telescope-seeding-config.json",
    )
    addCKFTracks(
        s,
        trackingGeometry,
        field,
        TrackSelectorConfig(),
        CkfConfig(
            chi2CutOffMeasurement=150.0,
            chi2CutOffOutlier=250.0,
            numMeasurementsCutOff=50,
            seedDeduplication=True,
            stayOnSeed=True,
        ),
        outputDirRoot=outputDir,
        writeTrackStates=True,
        writePerformance=False,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", "-o", type=Path, default=Path.cwd() / "misaligned_output"
    )
    parser.add_argument("--events", "-n", type=int, default=1000)
    parser.add_argument(
        "--layer",
        type=int,
        default=4,
        help="geo id layer to misalign (2, 4, ..., 18)",
    )
    parser.add_argument(
        "--shift",
        type=float,
        default=0.2,
        help="global z shift of that layer in mm",
    )
    parser.add_argument(
        "--target",
        choices=["sim", "reco", "both"],
        default="sim",
        help=(
            "which geometry context sees the shift. 'sim' is the interesting one: "
            "the detector is built shifted and reconstructed as if it were not. "
            "'both' means the shift is perfectly known and residuals stay centred."
        ),
    )
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    detector = TelescopeDetector(
        bounds=[200, 200],
        positions=LAYER_POSITIONS,
        stereos=[0] * len(LAYER_POSITIONS),
        binValue=1,
    )
    trackingGeometry = detector.trackingGeometry()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    s = Sequencer(events=args.events, numThreads=1, outputDir=str(args.output))

    addMisalignment(
        s,
        trackingGeometry,
        layer=args.layer,
        shift=acts.Vector3(0, 0, args.shift * u.mm),
        target={
            "sim": AlignmentDecorator.Target.eSim,
            "reco": AlignmentDecorator.Target.eReco,
            "both": AlignmentDecorator.Target.eBoth,
        }[args.target],
    )

    addTelescopeChain(s, trackingGeometry, field, args.output, RandomNumbers(seed=42))

    s.run()

    print(
        f"\nWrote {args.output / 'trackstates_ckf.root'}. The residual "
        f"res_eLOC0/res_eLOC1 for volume 1 layer {args.layer} carries the "
        f"{args.shift} mm shift when --target sim, and is centred on zero when "
        f"--target both."
    )


if __name__ == "__main__":
    main()
