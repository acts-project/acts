#!/usr/bin/env python3
from pathlib import Path
import argparse
import tempfile
import shutil
import os
import sys
import subprocess

sys.path += [
    str(Path(__file__).parent.parent.parent / "Examples/Scripts/Python/"),
]

# this has to happen before we import the ACTS module
import acts.examples

# @TODO: Fix failure in gain matrix smoothing
# See https://github.com/acts-project/acts/issues/1215
acts.logging.setFailureThreshold(acts.logging.FATAL)

from truth_tracking_kalman import runTruthTrackingKalman
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import (
    addParticleGun,
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
    SeedfinderConfigArg,
    SeedingAlgorithm,
    TrackParamsEstimationConfig,
    addCKFTracks,
    CKFPerformanceConfig,
)


parser = argparse.ArgumentParser()
parser.add_argument("outdir")
parser.add_argument("--events", type=int, default=10000)
parser.add_argument("--skip", type=int, default=0)

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

s = acts.examples.Sequencer(
    events=args.events, numThreads=-1, logLevel=acts.logging.INFO, skip=0
)

with tempfile.TemporaryDirectory() as temp:
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


### CKF track finding variations

for truthSmearedSeeded, truthEstimatedSeeded, label in [
    (True, False, "truth_smeared"),  # if first is true, second is ignored
    (False, True, "truth_estimated"),
    (False, False, "seeded"),
]:
    s = acts.examples.Sequencer(
        events=args.events, numThreads=1, logLevel=acts.logging.INFO, skip=args.skip
    )

    with tempfile.TemporaryDirectory() as temp:
        tp = Path(temp)

        for d in decorators:
            s.addContextDecorator(d)

        rnd = acts.examples.RandomNumbers(seed=42)

        s = addParticleGun(
            s,
            EtaConfig(-4.0, 4.0),
            ParticleConfig(4, acts.PdgParticle.eMuon, True),
            PhiConfig(0.0, 360.0 * u.degree),
            multiplicity=2,
            rnd=rnd,
        )

        s = addFatras(
            s,
            trackingGeometry,
            field,
            rnd=rnd,
        )

        s = addDigitization(
            s,
            trackingGeometry,
            field,
            digiConfigFile=digiConfig,
            rnd=rnd,
        )

        s = addSeeding(
            s,
            trackingGeometry,
            field,
            TruthSeedRanges(pt=(500.0 * u.MeV, None), nHits=(9, None)),
            ParticleSmearingSigmas(
                pRel=0.01
            ),  # only used by SeedingAlgorithm.TruthSmeared
            SeedfinderConfigArg(
                r=(None, 200 * u.mm),  # rMin=default, 33mm
                deltaR=(1 * u.mm, 60 * u.mm),
                collisionRegion=(-250 * u.mm, 250 * u.mm),
                z=(-2000 * u.mm, 2000 * u.mm),
                maxSeedsPerSpM=1,
                sigmaScattering=50,
                radLengthPerSeed=0.1,
                minPt=500 * u.MeV,
                bFieldInZ=1.99724 * u.T,
                impactMax=3 * u.mm,
            ),
            TrackParamsEstimationConfig(deltaR=(10.0 * u.mm, None)),
            seedingAlgorithm=SeedingAlgorithm.TruthSmeared
            if truthSmearedSeeded
            else SeedingAlgorithm.TruthEstimated
            if truthEstimatedSeeded
            else SeedingAlgorithm.Default,
            geoSelectionConfigFile=geoSel,
            outputDirRoot=tp,
            rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
        )

        s = addCKFTracks(
            s,
            trackingGeometry,
            field,
            CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
            outputDirRoot=tp,
            outputDirCsv=None,
        )

        s.run()
        del s

        perf_file = tp / "performance_ckf.root"
        assert perf_file.exists(), "Performance file not found"
        shutil.copy(perf_file, outdir / f"performance_ckf_tracks_{label}.root")

        if not truthSmearedSeeded and not truthEstimatedSeeded:
            residual_app = srcdir / "build/bin/ActsAnalysisResidualsAndPulls"
            # @TODO: Add try/except
            subprocess.check_call(
                [
                    str(residual_app),
                    "--predicted",
                    "--filtered",
                    "--smoothed",
                    "--silent",
                    "-i",
                    str(tp / "trackstates_ckf.root"),
                    "-o",
                    str(outdir / "acts_analysis_residuals_and_pulls.root"),
                    "--save",
                    "",
                ]
            )
