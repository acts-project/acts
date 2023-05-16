#!/usr/bin/env python3
from pathlib import Path
import argparse
import tempfile
import shutil
import datetime
import sys

sys.path += [
    str(Path(__file__).parent.parent.parent / "Examples/Scripts/Python/"),
]

# this has to happen before we import the ACTS module
import acts.examples

# @TODO: Fix failure in gain matrix smoothing
# See https://github.com/acts-project/acts/issues/1215
acts.logging.setFailureThreshold(acts.logging.FATAL)

from truth_tracking_kalman import runTruthTrackingKalman
from truth_tracking_gsf import runTruthTrackingGsf
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    ParticleSmearingSigmas,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    TruthEstimatedSeedingAlgorithmConfigArg,
    addCKFTracks,
    CKFPerformanceConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    TrackSelectorRanges,
)


parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=["all", "kalman", "gsf", "fullchains", "vertexing"])
parser.add_argument("outdir")

args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(exist_ok=True)


srcdir = Path(__file__).resolve().parent.parent.parent


u = acts.UnitConstants

matDeco = acts.IMaterialDecorator.fromFile(
    srcdir / "thirdparty/OpenDataDetector/data/odd-material-maps.root",
    level=acts.logging.INFO,
)
detector, trackingGeometry, decorators = getOpenDataDetector(
    getOpenDataDetectorDirectory(), matDeco
)
digiConfig = srcdir / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json"
geoSel = srcdir / "thirdparty/OpenDataDetector/config/odd-seeding-config.json"

field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


def truth_tracking_kalman():
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=10000, numThreads=-1, logLevel=acts.logging.INFO
        )

        tp = Path(temp)
        runTruthTrackingKalman(
            trackingGeometry,
            field,
            digiConfigFile=digiConfig,
            outputDir=tp,
            s=s,
        )

        s.run()
        del s

        perf_file = tp / "performance_track_fitter.root"
        assert perf_file.exists(), "Performance file not found"
        shutil.copy(perf_file, outdir / "performance_truth_tracking.root")


def truth_tracking_gsf():
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=500, numThreads=-1, logLevel=acts.logging.INFO
        )

        tp = Path(temp)
        runTruthTrackingGsf(
            trackingGeometry,
            digiConfig,
            field,
            outputDir=tp,
            s=s,
        )

        s.run()
        del s

        perf_file = tp / "performance_gsf.root"
        assert perf_file.exists(), "Performance file not found"
        shutil.copy(perf_file, outdir / "performance_gsf.root")


def run_ckf_tracking(truthSmearedSeeded, truthEstimatedSeeded, label):
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=500, numThreads=-1, logLevel=acts.logging.INFO
        )

        tp = Path(temp)

        for d in decorators:
            s.addContextDecorator(d)

        rnd = acts.examples.RandomNumbers(seed=42)

        addParticleGun(
            s,
            EtaConfig(-4.0, 4.0),
            ParticleConfig(4, acts.PdgParticle.eMuon, True),
            PhiConfig(0.0, 360.0 * u.degree),
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(10 * u.um, 10 * u.um, 50 * u.mm, 0),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=50,
            rnd=rnd,
        )

        addFatras(
            s,
            trackingGeometry,
            field,
            rnd=rnd,
        )

        addDigitization(
            s,
            trackingGeometry,
            field,
            digiConfigFile=digiConfig,
            rnd=rnd,
        )

        addSeeding(
            s,
            trackingGeometry,
            field,
            TruthSeedRanges(pt=(500.0 * u.MeV, None), nHits=(9, None)),
            ParticleSmearingSigmas(
                pRel=0.01
            ),  # only used by SeedingAlgorithm.TruthSmeared
            SeedFinderConfigArg(
                r=(None, 200 * u.mm),  # rMin=default, 33mm
                deltaR=(1 * u.mm, 60 * u.mm),
                collisionRegion=(-250 * u.mm, 250 * u.mm),
                z=(-2000 * u.mm, 2000 * u.mm),
                maxSeedsPerSpM=1,
                sigmaScattering=5,
                radLengthPerSeed=0.1,
                minPt=500 * u.MeV,
                impactMax=3 * u.mm,
            ),
            SeedFinderOptionsArg(bFieldInZ=1.99724 * u.T, beamPos=(0.0, 0.0)),
            TruthEstimatedSeedingAlgorithmConfigArg(deltaR=(10.0 * u.mm, None)),
            seedingAlgorithm=SeedingAlgorithm.TruthSmeared
            if truthSmearedSeeded
            else SeedingAlgorithm.TruthEstimated
            if truthEstimatedSeeded
            else SeedingAlgorithm.Default
            if label == "seeded"
            else SeedingAlgorithm.Orthogonal,
            geoSelectionConfigFile=geoSel,
            rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
            outputDirRoot=tp,
        )

        addCKFTracks(
            s,
            trackingGeometry,
            field,
            CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
            TrackSelectorRanges(
                loc0=(-4.0 * u.mm, 4.0 * u.mm),
                pt=(500 * u.MeV, None),
            ),
            outputDirRoot=tp,
        )

        if label in ["seeded", "orthogonal"]:
            addAmbiguityResolution(
                s,
                AmbiguityResolutionConfig(maximumSharedHits=3),
                CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
                outputDirRoot=tp,
            )

        addVertexFitting(
            s,
            field,
            associatedParticles=None
            if label in ["seeded", "orthogonal"]
            else "particles_input",
            outputProtoVertices="ivf_protovertices",
            outputVertices="ivf_fittedVertices",
            vertexFinder=VertexFinder.Iterative,
            outputDirRoot=tp / "ivf",
        )

        addVertexFitting(
            s,
            field,
            associatedParticles=None
            if label in ["seeded", "orthogonal"]
            else "particles_input",
            outputProtoVertices="amvf_protovertices",
            outputVertices="amvf_fittedVertices",
            outputTime="amvf_outputTime",
            vertexFinder=VertexFinder.AMVF,
            outputDirRoot=tp / "amvf",
        )

        s.run()
        del s

        for vertexing in ["ivf", "amvf"]:
            shutil.move(
                tp / f"{vertexing}/performance_vertexing.root",
                tp / f"performance_{vertexing}.root",
            )

        for stem in [
            "performance_ckf",
            "tracksummary_ckf",
            "performance_ivf",
            "performance_amvf",
        ] + (
            ["performance_seeding", "performance_ambi"]
            if label in ["seeded", "orthogonal"]
            else ["performance_seeding"]
            if label == "truth_estimated"
            else []
        ):
            perf_file = tp / f"{stem}.root"
            assert perf_file.exists(), "Performance file not found"
            shutil.copy(perf_file, outdir / f"{stem}_{label}.root")


def run_vertexing(fitter, mu, events):
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=events, numThreads=-1, logLevel=acts.logging.INFO
        )

        tp = Path(temp)

        for d in decorators:
            s.addContextDecorator(d)

        rnd = acts.examples.RandomNumbers(seed=42)

        addParticleGun(
            s,
            MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
            EtaConfig(-3.0, 3.0),
            ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(10 * u.um, 10 * u.um, 50 * u.mm, 0),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=mu,
            rnd=rnd,
        )

        addFatras(
            s,
            trackingGeometry,
            field,
            rnd=rnd,
        )

        addDigitization(
            s,
            trackingGeometry,
            field,
            digiConfigFile=digiConfig,
            rnd=rnd,
        )

        addSeeding(
            s,
            trackingGeometry,
            field,
            SeedFinderConfigArg(
                r=(None, 200 * u.mm),  # rMin=default, 33mm
                deltaR=(1 * u.mm, 60 * u.mm),
                collisionRegion=(-250 * u.mm, 250 * u.mm),
                z=(-2000 * u.mm, 2000 * u.mm),
                maxSeedsPerSpM=1,
                sigmaScattering=5,
                radLengthPerSeed=0.1,
                minPt=500 * u.MeV,
                impactMax=3 * u.mm,
            ),
            SeedFinderOptionsArg(bFieldInZ=1.99724 * u.T),
            seedingAlgorithm=SeedingAlgorithm.Default,
            geoSelectionConfigFile=geoSel,
            rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
        )

        addCKFTracks(
            s,
            trackingGeometry,
            field,
            CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
            TrackSelectorRanges(
                pt=(500 * u.MeV, None),
                loc0=(-4.0 * u.mm, 4.0 * u.mm),
                absEta=(None, 3.0),
            ),
        )

        addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(maximumSharedHits=3),
            CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
        )

        addVertexFitting(
            s,
            field,
            vertexFinder=fitter,
            outputDirRoot=tp,
        )

        s.run()

        del s

        perf_file = tp / f"performance_vertexing.root"
        assert perf_file.exists(), "Performance file not found"
        shutil.copy(
            perf_file,
            outdir / f"performance_vertexing_{fitter.name}_mu{mu}.root",
        )


with acts.FpeMonitor():

    ### Truth tracking with Kalman Filter
    if args.mode == "all" or args.mode == "kalman":
        truth_tracking_kalman()

    ### GSF
    if args.mode == "all" or args.mode == "gsf":
        truth_tracking_gsf()

    ### CKF track finding variations
    if args.mode == "all" or args.mode == "fullchains":
        for truthSmearedSeeded, truthEstimatedSeeded, label in [
            (True, False, "truth_smeared"),  # if first is true, second is ignored
            (False, True, "truth_estimated"),
            (False, False, "seeded"),
            (False, False, "orthogonal"),
        ]:
            run_ckf_tracking(truthSmearedSeeded, truthEstimatedSeeded, label)

    ### VERTEX MU SCAN
    if args.mode == "all" or args.mode == "vertexing":
        for fitter in (VertexFinder.Iterative, VertexFinder.AMVF):
            for mu in (1, 10, 25, 50, 75, 100, 125, 150, 175, 200):
                start = datetime.datetime.now()

                events = 5
                run_vertexing(fitter, mu, events)

                delta = datetime.datetime.now() - start

                duration = delta.total_seconds() / events

                (
                    outdir / f"performance_vertexing_{fitter.name}_mu{mu}_time.txt"
                ).write_text(str(duration))
