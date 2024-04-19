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
    CkfConfig,
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
            events=500,
            numThreads=-1,
            logLevel=acts.logging.INFO,
        )

        tp = Path(temp)

        for d in setup.decorators:
            s.addContextDecorator(d)

        rnd = acts.examples.RandomNumbers(seed=42)

        addParticleGun(
            s,
            MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
            EtaConfig(-3.0, 3.0),
            PhiConfig(0.0, 360.0 * u.degree),
            ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns
                ),
            ),
            multiplicity=50,
            rnd=rnd,
        )

        addFatras(
            s,
            setup.trackingGeometry,
            setup.field,
            enableInteractions=True,
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
            TruthSeedRanges(pt=(500 * u.MeV, None), nHits=(9, None)),
            ParticleSmearingSigmas(
                pRel=0.01
            ),  # only used by SeedingAlgorithm.TruthSmeared
            SeedFinderConfigArg(
                r=(33 * u.mm, 200 * u.mm),
                deltaR=(1 * u.mm, 60 * u.mm),
                collisionRegion=(-250 * u.mm, 250 * u.mm),
                z=(-2000 * u.mm, 2000 * u.mm),
                maxSeedsPerSpM=1,
                sigmaScattering=5,
                radLengthPerSeed=0.1,
                minPt=500 * u.MeV,
                impactMax=3 * u.mm,
            ),
            SeedFinderOptionsArg(bFieldInZ=2 * u.T),
            TruthEstimatedSeedingAlgorithmConfigArg(deltaR=(10.0 * u.mm, None)),
            seedingAlgorithm=SeedingAlgorithm.TruthSmeared
            if truthSmearedSeeded
            else SeedingAlgorithm.TruthEstimated
            if truthEstimatedSeeded
            else SeedingAlgorithm.Default
            if label == "seeded"
            else SeedingAlgorithm.Orthogonal,
            initialSigmas=[
                1 * u.mm,
                1 * u.mm,
                1 * u.degree,
                1 * u.degree,
                0.1 / u.GeV,
                1 * u.ns,
            ],
            initialVarInflation=[1.0] * 6,
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
                maxHoles=2,
                maxOutliers=2,
            ),
            CkfConfig(
                seedDeduplication=False if truthSmearedSeeded else True,
                stayOnSeed=False if truthSmearedSeeded else True,
            ),
            outputDirRoot=tp,
        )

        if label in ["seeded", "orthogonal"]:
            addAmbiguityResolution(
                s,
                AmbiguityResolutionConfig(
                    maximumSharedHits=3,
                    maximumIterations=10000,
                    nMeasurementsMin=6,
                ),
                outputDirRoot=tp,
            )

        s.addAlgorithm(
            acts.examples.TracksToParameters(
                level=acts.logging.INFO,
                inputTracks="tracks",
                outputTrackParameters="trackParameters",
            )
        )

        # Choosing a seeder only has an effect on VertexFinder.AMVF. For
        # VertexFinder.IVF we always use acts.VertexSeedFinder.GaussianSeeder
        # (Python binding is not implemented).
        # Setting useTime also only has an effect on VertexFinder.AMVF due to
        # the same reason.
        addVertexFitting(
            s,
            setup.field,
            trackParameters="trackParameters",
            outputProtoVertices="ivf_protovertices",
            outputVertices="ivf_fittedVertices",
            vertexFinder=VertexFinder.Iterative,
            outputDirRoot=tp / "ivf",
        )

        addVertexFitting(
            s,
            setup.field,
            trackParameters="trackParameters",
            outputProtoVertices="amvf_protovertices",
            outputVertices="amvf_fittedVertices",
            seeder=acts.VertexSeedFinder.GaussianSeeder,
            useTime=False,  # Time seeding not implemented for the Gaussian seeder
            vertexFinder=VertexFinder.AMVF,
            outputDirRoot=tp / "amvf",
        )

        # Use the adaptive grid vertex seeder in combination with the AMVF
        # To avoid having too many physmon cases, we only do this for the label
        # "seeded"
        if label == "seeded":
            addVertexFitting(
                s,
                setup.field,
                trackParameters="trackParameters",
                outputProtoVertices="amvf_gridseeder_protovertices",
                outputVertices="amvf_gridseeder_fittedVertices",
                seeder=acts.VertexSeedFinder.AdaptiveGridSeeder,
                useTime=True,
                vertexFinder=VertexFinder.AMVF,
                outputDirRoot=tp / "amvf_gridseeder",
            )

        s.run()
        del s

        for vertexing in ["ivf", "amvf"]:
            shutil.move(
                tp / f"{vertexing}/performance_vertexing.root",
                tp / f"performance_{vertexing}.root",
            )

        if label == "seeded":
            vertexing = "amvf_gridseeder"
            shutil.move(
                tp / f"{vertexing}/performance_vertexing.root",
                tp / f"performance_{vertexing}.root",
            )

        for stem in (
            [
                "performance_ckf",
                "tracksummary_ckf",
                "performance_ivf",
                "performance_amvf",
            ]
            + (["performance_amvf_gridseeder"] if label == "seeded" else [])
            + (
                ["performance_seeding", "performance_ambi"]
                if label in ["seeded", "orthogonal"]
                else ["performance_seeding"]
                if label == "truth_estimated"
                else []
            )
        ):
            perf_file = tp / f"{stem}.root"
            assert perf_file.exists(), "Performance file not found"
            shutil.copy(perf_file, setup.outdir / f"{stem}_{label}.root")


for truthSmearedSeeded, truthEstimatedSeeded, label in [
    (True, False, "truth_smeared"),  # if first is true, second is ignored
    (False, True, "truth_estimated"),
    (False, False, "seeded"),
    (False, False, "orthogonal"),
]:
    run_ckf_tracking(truthSmearedSeeded, truthEstimatedSeeded, label)
