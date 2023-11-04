#!/usr/bin/env python3
import tempfile
from pathlib import Path
import shutil

import acts
from acts.examples.simulation import (
    addPythia8,
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


with tempfile.TemporaryDirectory() as temp:
    fpeMasks = acts.examples.Sequencer.FpeMask.fromFile(
        Path(__file__).parent.parent / "fpe_masks.yml"
    ) + [
        acts.examples.Sequencer.FpeMask(
            "Examples/Algorithms/Fatras/src/FatrasSimulation.cpp",
            (172, 173),
            acts.FpeType.FLTINV,
            1,
        ),
        acts.examples.Sequencer.FpeMask(
            "Examples/Algorithms/Fatras/src/FatrasSimulation.cpp",
            (235, 236),
            acts.FpeType.FLTINV,
            1,
        ),
        acts.examples.Sequencer.FpeMask(
            "Examples/Algorithms/Fatras/src/FatrasSimulation.cpp",
            (235, 236),
            acts.FpeType.FLTOVF,
            1,
        ),
        acts.examples.Sequencer.FpeMask(
            "Examples/Io/Root/src/RootTrajectorySummaryWriter.cpp",
            (371, 372),
            acts.FpeType.FLTINV,
            1,
        ),
        acts.examples.Sequencer.FpeMask(
            "Core/src/Utilities/AnnealingUtility.cpp",
            (38, 39),
            acts.FpeType.FLTUND,
            1,
        ),
        acts.examples.Sequencer.FpeMask(
            "Fatras/include/ActsFatras/Kernel/detail/SimulationActor.hpp",
            (110, 111),
            acts.FpeType.FLTINV,
            1,
        ),
        acts.examples.Sequencer.FpeMask(
            "Fatras/include/ActsFatras/Kernel/Simulation.hpp",
            (96, 97),
            acts.FpeType.FLTOVF,
            1,
        ),
    ]
    s = acts.examples.Sequencer(
        events=3,
        numThreads=-1,
        logLevel=acts.logging.INFO,
        fpeMasks=fpeMasks,
    )

    tp = Path(temp)

    for d in setup.decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    addPythia8(
        s,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=200,
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
        ),
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
        ParticleSmearingSigmas(pRel=0.01),  # only used by SeedingAlgorithm.TruthSmeared
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
        seedingAlgorithm=SeedingAlgorithm.Default,
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

    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(maximumSharedHits=3),
        outputDirRoot=tp,
    )

    addVertexFitting(
        s,
        setup.field,
        seeder=acts.VertexSeedFinder.GaussianSeeder,
        outputProtoVertices="amvf_protovertices",
        outputVertices="amvf_fittedVertices",
        vertexFinder=VertexFinder.AMVF,
        outputDirRoot=tp / "amvf",
    )

    addVertexFitting(
        s,
        setup.field,
        seeder=acts.VertexSeedFinder.AdaptiveGridSeeder,
        outputProtoVertices="amvf_gridseeder_protovertices",
        outputVertices="amvf_gridseeder_fittedVertices",
        vertexFinder=VertexFinder.AMVF,
        outputDirRoot=tp / "amvf_gridseeder",
    )

    s.run()
    del s

    for vertexing in ["amvf", "amvf_gridseeder"]:
        shutil.move(
            tp / f"{vertexing}/performance_vertexing.root",
            tp / f"performance_{vertexing}.root",
        )

    for stem in [
        "performance_ckf",
        "tracksummary_ckf",
        "performance_amvf",
        "performance_amvf_gridseeder",
    ] + (["performance_seeding", "performance_ambi"]):
        perf_file = tp / f"{stem}.root"
        assert perf_file.exists(), "Performance file not found"
        shutil.copy(perf_file, setup.outdir / f"{stem}_ttbar.root")
