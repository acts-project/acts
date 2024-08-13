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
    addTracccChain,
)

from physmon_common import makeSetup

setup = makeSetup(useGeometricConfig=True)


def runTracccChain(
    platform,
    actsDigi=False,
    truthSmearedSeeded=False,
    truthEstimatedSeeded=False,
    label="seeded"):

    u = acts.UnitConstants
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=500,
            numThreads=-1,
            logLevel=acts.logging.INFO,
            trackFpes=False,
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
            minEnergyDeposit=0,
            # Traccc chain expects doMerge and mergeCommonCorner 
            # as ACTS and Traccc clusterization must give the 
            # same results to be able to pair the measurements
            # (which is used for truth matching).
            doMerge=True,
            mergeCommonCorner=True,
        )

        if actsDigi:
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
                seedingAlgorithm=(
                    SeedingAlgorithm.TruthSmeared
                    if truthSmearedSeeded
                    else (
                        SeedingAlgorithm.TruthEstimated
                        if truthEstimatedSeeded
                        else (
                            SeedingAlgorithm.Default
                            if label == "seeded"
                            else SeedingAlgorithm.Orthogonal
                        )
                    )
                ),
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
            )

        chainConfig = acts.examples.TracccChainConfig()

        addTracccChain(
            s,
            setup.trackingGeometry,
            setup.field,
            digiConfigFile=setup.digiConfig,
            inputCells="cells",
            outputDirRoot=tp,
            chainConfig=chainConfig,
            logLevel=acts.logging.INFO,
            platform=platform,
            enableAmbiguityResolution=True,
            reconstructionOnly=actsDigi,
        )

        s.run()
        del s

        recStr = "_reconstruction_only" if actsDigi else ""
        name = f"tracksummary_traccc_{platform + recStr}"
        perf_file = tp / f"{name}.root"
        assert perf_file.exists(), "Tracksummary file not found"
        shutil.copy(perf_file, setup.outdir / f"{name}.root")

        if not actsDigi:
            name = f"performance_seeding_traccc_{platform}";
            perf_file = tp / f"{name}.root"
            assert perf_file.exists(), "Performance seeding file not found"
            shutil.copy(perf_file, setup.outdir / f"{name}.root")