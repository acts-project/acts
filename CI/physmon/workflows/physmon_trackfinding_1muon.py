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
    TrackSmearingSigmas,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    GbtsConfigArg,
    CkfConfig,
    addCKFTracks,
    TrackSelectorConfig,
)

from physmon_common import (
    copyOutputs,
    GBTS_CONNECTION_TABLE,
    GBTS_LAYER_MAP,
    SEED_PARAMETER_SIGMAS,
    makeSetup,
    seedingVariants,
)

u = acts.UnitConstants

# Seeding variants, as `SeedingAlgorithm` member names (see seedingVariants).
# config/physmon.yaml picks which of them run in CI.
VARIANTS = [
    "TruthSmeared",
    "TruthEstimated",
    "GridTriplet",
    "OrthogonalTriplet",
    "Gbts",
]

setup = makeSetup(VARIANTS)

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
    seedFinderOptionsArg=SeedFinderOptionsArg(bFieldInZ=setup.bFieldInZ),
    geoSelectionConfigFile=setup.geoSel,
    **SEED_PARAMETER_SIGMAS,
)

CKF = dict(
    chi2CutOffMeasurement=15.0,
    chi2CutOffOutlier=25.0,
    numMeasurementsCutOff=10,
)

# The files the chain writes, and the names they are stored under.
OUTPUTS = {
    "performance_seeding.root": "performance_seeding.root",
    "tracksummary_ckf.root": "tracksummary_ckf.root",
    "performance_finding_ckf.root": "performance_finding_ckf.root",
    "performance_fitting_ckf.root": "performance_fitting_ckf.root",
}


def seedingSettings(seeding):
    """The addSeeding arguments that only this seeding algorithm needs."""
    spreadTriplet = acts.examples.SeedSpacePointSelection.SpreadTriplet
    if seeding == SeedingAlgorithm.TruthSmeared:
        return dict(
            # zero everything so the CKF has a chance to find the measurements
            trackSmearingSigmas=TrackSmearingSigmas(
                loc0=0,
                loc0PtA=0,
                loc0PtB=0,
                loc1=0,
                loc1PtA=0,
                loc1PtB=0,
                time=0,
                phi=0,
                theta=0,
                ptRel=0,
            ),
        )
    if seeding == SeedingAlgorithm.TruthEstimated:
        return dict(paramEstimationSpacePoints=spreadTriplet)
    if seeding == SeedingAlgorithm.Gbts:
        return dict(
            paramEstimationSpacePoints=spreadTriplet,
            layerMappingConfigFile=GBTS_LAYER_MAP,
            connectorInputConfigFile=GBTS_CONNECTION_TABLE,
            gbtsConfigArg=GbtsConfigArg(
                seederConfig=acts.examples.GraphBasedSeedingConfig(
                    minSeedLevel=2,
                ),
                graphConfig=acts.examples.GbtsGraphBuilderConfig(
                    minPt=0.9 * u.GeV,
                    z0HistogramMaxBarrelOrder=-1,
                ),
                useAdaptiveCuts=False,
            ),
            gbtsSeedParameters=True,
        )
    return {}


def run_ckf_tracking(seeding, outdir):
    # TruthSmeared starts the CKF from smeared truth parameters, without seeds
    hasSeeds = seeding != SeedingAlgorithm.TruthSmeared

    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=10000,
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
            ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns
                ),
            ),
            multiplicity=1,
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
            seedingAlgorithm=seeding,
            rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
            outputDirRoot=tp,
            **(SEEDING | seedingSettings(seeding)),
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
            CkfConfig(**CKF, seedDeduplication=hasSeeds, stayOnSeed=hasSeeds),
            outputDirRoot=tp,
        )

        s.run()

        outputs = dict(OUTPUTS)
        if not hasSeeds:
            del outputs["performance_seeding.root"]
        copyOutputs(outputs, tp, outdir)


for seeding, outdir in seedingVariants(setup, VARIANTS):
    run_ckf_tracking(seeding, outdir)
