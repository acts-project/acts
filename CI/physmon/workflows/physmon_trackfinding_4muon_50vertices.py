#!/usr/bin/env python3

import tempfile
from pathlib import Path

import acts
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
    ParticleSelectorConfig,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    CkfConfig,
    addCKFTracks,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    TrackSelectorConfig,
)

from physmon_common import SEED_PARAMETER_SIGMAS, copyOutputs, makeSetup

u = acts.UnitConstants

setup = makeSetup()

SEEDING = dict(
    seedFinderConfigArg=SeedFinderConfigArg(
        r=(33 * u.mm, 200 * u.mm),
        deltaR=(1 * u.mm, 300 * u.mm),
        collisionRegion=(-250 * u.mm, 250 * u.mm),
        z=(-2000 * u.mm, 2000 * u.mm),
        maxSeedsPerSpM=1,
        sigmaScattering=5,
        radLengthPerSeed=0.1,
        minPt=0.5 * u.GeV,
        impactMax=3 * u.mm,
    ),
    seedFinderOptionsArg=SeedFinderOptionsArg(
        bFieldInZ=setup.bFieldInZ, beamPos=(0.0, 0.0)
    ),
    geoSelectionConfigFile=setup.geoSel,
    **SEED_PARAMETER_SIGMAS,
)

CKF = dict(
    chi2CutOffMeasurement=15.0,
    chi2CutOffOutlier=25.0,
    numMeasurementsCutOff=10,
    seedDeduplication=True,
    stayOnSeed=True,
)

# Vertex finders, keyed by the name of their output directory and file.
# The seeder and useTime only have an effect on VertexFinder.AMVF: VertexFinder.IVF
# always uses acts.examples.VertexSeedFinder.GaussianSeeder (Python binding is
# not implemented).
VERTEXING = {
    "ivf_notime": dict(vertexFinder=VertexFinder.Iterative),
    "amvf_gauss_notime": dict(
        vertexFinder=VertexFinder.AMVF,
        seeder=acts.examples.VertexSeedFinder.GaussianSeeder,
        useTime=False,  # Time seeding not implemented for the Gaussian seeder
    ),
    "amvf_grid_time": dict(
        vertexFinder=VertexFinder.AMVF,
        seeder=acts.examples.VertexSeedFinder.AdaptiveGridSeeder,
        useTime=True,
    ),
}

# The files the chain writes, and the names they are stored under.
OUTPUTS = {
    "performance_seeding.root": "performance_seeding.root",
    "tracksummary_ckf.root": "tracksummary_ckf.root",
    "performance_finding_ckf.root": "performance_finding_ckf.root",
    "performance_fitting_ckf.root": "performance_fitting_ckf.root",
    "performance_finding_ambi.root": "performance_finding_ckf_ambi.root",
    "performance_fitting_ambi.root": "performance_fitting_ckf_ambi.root",
    "ivf_notime/performance_vertexing.root": "performance_vertexing_ivf_notime.root",
    "amvf_gauss_notime/performance_vertexing.root": "performance_vertexing_amvf_gauss_notime.root",
    "amvf_grid_time/performance_vertexing.root": "performance_vertexing_amvf_grid_time.root",
}

with tempfile.TemporaryDirectory() as temp:
    s = acts.examples.Sequencer(
        events=3,
        numThreads=setup.threads,
        logLevel=acts.logging.INFO,
    )

    tp = Path(temp)

    for d in setup.decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-3.0, 3.0, uniform=True),
        PhiConfig(0.0, 360.0 * u.degree),
        ParticleConfig(10, acts.PdgParticle.eMuon, randomizeCharge=True),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns),
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

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(0.9 * u.GeV, None),
            measurements=(9, None),
            removeNeutral=True,
        ),
    )

    addSeeding(
        s,
        setup.trackingGeometry,
        setup.field,
        seedingAlgorithm=SeedingAlgorithm.GridTriplet,
        outputDirRoot=tp,
        **SEEDING,
    )

    addCKFTracks(
        s,
        setup.trackingGeometry,
        setup.field,
        TrackSelectorConfig(
            pt=(0.9 * u.GeV, None),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
            nMeasurementsMin=6,
            maxHoles=2,
            maxOutliers=2,
        ),
        CkfConfig(**CKF),
        outputDirRoot=tp,
    )

    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
            maximumSharedHits=3,
            maximumIterations=100000,
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

    for name, settings in VERTEXING.items():
        addVertexFitting(
            s,
            setup.field,
            trackParameters="trackParameters",
            outputProtoVertices=f"{name}_protovertices",
            outputVertices=f"{name}_vertices",
            outputVertexTruthMatching=f"{name}_vertex_truth_matching",
            outputTruthVertexMatching=f"{name}_truth_vertex_matching",
            outputDirRoot=tp / name,
            writeTrackInfo=True,
            **settings,
        )

    s.run()

    copyOutputs(OUTPUTS, tp, setup.outdir)
