#!/usr/bin/env python3
import tempfile
from pathlib import Path
import shutil
import datetime

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
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
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


def run_vertexing(fitter, mu, events):
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=events,
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
            multiplicity=mu,
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
            seedingAlgorithm=SeedingAlgorithm.Default,
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
        )

        addCKFTracks(
            s,
            setup.trackingGeometry,
            setup.field,
            TrackSelectorConfig(
                loc0=(-4.0 * u.mm, 4.0 * u.mm),
                pt=(500 * u.MeV, None),
                nMeasurementsMin=6,
                maxHoles=2,
                maxOutliers=2,
            ),
            CkfConfig(
                seedDeduplication=True,
                stayOnSeed=True,
            ),
        )

        addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(
                maximumSharedHits=3,
                maximumIterations=10000,
                nMeasurementsMin=6,
            ),
        )

        addVertexFitting(
            s,
            setup.field,
            vertexFinder=fitter,
            outputDirRoot=tp,
        )

        s.run()

        del s

        perf_file = tp / f"performance_vertexing.root"
        assert perf_file.exists(), "Performance file not found"
        shutil.copy(
            perf_file,
            setup.outdir / f"performance_vertexing_{fitter.name}_mu{mu}.root",
        )


for fitter in (VertexFinder.Iterative, VertexFinder.AMVF):
    for mu in (1, 10, 25, 50, 75, 100, 125, 150, 175, 200):
        start = datetime.datetime.now()

        events = 5
        run_vertexing(fitter, mu, events)

        delta = datetime.datetime.now() - start

        duration = delta.total_seconds() / events

        (
            setup.outdir / f"performance_vertexing_{fitter.name}_mu{mu}_time.txt"
        ).write_text(str(duration))
