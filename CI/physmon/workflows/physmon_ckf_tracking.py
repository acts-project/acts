#!/usr/bin/env python3
import tempfile
from pathlib import Path
import shutil

import acts
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
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    TrackSelectorConfig,
)

from physmon_common import makeSetup

u = acts.UnitConstants

setup = makeSetup()


def run_ckf_tracking(truthSmearedSeeded, truthEstimatedSeeded, label):
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=500, numThreads=-1, logLevel=acts.logging.INFO
        )

        tp = Path(temp)

        for d in setup.decorators:
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
            setup.trackingGeometry,
            setup.field,
            rnd=rnd,
        )

        addDigitization(
            s,
            setup.trackingGeometry,
            setup.field,
            digiConfigFile=setup.digiConfig,
            rnd=rnd,
        )

        addSeeding(
            s,
            setup.trackingGeometry,
            setup.field,
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
            geoSelectionConfigFile=setup.geoSel,
            rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
            outputDirRoot=tp,
        )

        addCKFTracks(
            s,
            setup.trackingGeometry,
            setup.field,
            TrackSelectorConfig(
                pt=(500 * u.MeV, None),
                loc0=(-4.0 * u.mm, 4.0 * u.mm),
                nMeasurementsMin=6,
            ),
            outputDirRoot=tp,
        )

        if label in ["seeded", "orthogonal"]:
            addAmbiguityResolution(
                s,
                AmbiguityResolutionConfig(maximumSharedHits=3),
                outputDirRoot=tp,
            )

        addVertexFitting(
            s,
            setup.field,
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
            setup.field,
            associatedParticles=None
            if label in ["seeded", "orthogonal"]
            else "particles_input",
            outputProtoVertices="amvf_protovertices",
            outputVertices="amvf_fittedVertices",
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
            shutil.copy(perf_file, setup.outdir / f"{stem}_{label}.root")


with acts.FpeMonitor():
    for truthSmearedSeeded, truthEstimatedSeeded, label in [
        (True, False, "truth_smeared"),  # if first is true, second is ignored
        (False, True, "truth_estimated"),
        (False, False, "seeded"),
        (False, False, "orthogonal"),
    ]:
        run_ckf_tracking(truthSmearedSeeded, truthEstimatedSeeded, label)
