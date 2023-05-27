#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil
from enum import Enum
import itertools

import acts
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addGeant4,
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
)

from physmon_common import makeSetup

u = acts.UnitConstants

setup = makeSetup()

Particle = Enum(
    "Particle",
    [acts.PdgParticle.eMuon, acts.PdgParticle.ePionPlus, acts.PdgParticle.eElectron],
)
PT = Enum("PT", [1 * u.GeV, 10 * u.GeV, 100 * u.GeV])
Simulation = Enum("Simulation", ["Fatras", "Geant4"])
Seeding = Enum(
    "Simulation",
    [
        SeedingAlgorithm.TruthSmeared,
        SeedingAlgorithm.TruthEstimated,
        SeedingAlgorithm.Default,
        SeedingAlgorithm.Orthogonal,
    ],
)


def run_single_particles(particle, pT, simulation, seeding, label):
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=10000, numThreads=-1, logLevel=acts.logging.INFO
        )

        tp = Path(temp)

        for d in setup.decorators:
            s.addContextDecorator(d)

        rnd = acts.examples.RandomNumbers(seed=42)

        addParticleGun(
            s,
            MomentumConfig(pT, pT, transverse=True),
            EtaConfig(-3.0, 3.0),
            ParticleConfig(1, particle, True),
            PhiConfig(0.0 * u.degree, 360.0 * u.degree),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(10 * u.um, 10 * u.um, 50 * u.mm, 1 * u.ns),
            ),
            multiplicity=1,
            rnd=rnd,
        )

        if simulation == Simulation.Fatras:
            addFatras(
                s,
                setup.trackingGeometry,
                setup.field,
                rnd=rnd,
            )
        elif simulation == Simulation.Geant4:
            addGeant4(
                s,
                setup.trackingGeometry,
                setup.field,
                rnd=rnd,
            )
        else:
            raise ValueError(f"unhandled simulation {simulation}")

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
            seedingAlgorithm=seeding,
            geoSelectionConfigFile=setup.geoSel,
            rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
            outputDirRoot=tp,
        )

        addCKFTracks(
            s,
            setup.trackingGeometry,
            setup.field,
            outputDirRoot=tp,
        )

        s.run()
        del s

        for stem in [
            "performance_ckf",
            "tracksummary_ckf",
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


def create_label(particle, pT, simulation, seeding):
    return f"{particle}_{pT}_{simulation}_{seeding}"


with acts.FpeMonitor():
    for particle, pT, simulation, seeding in itertools.product(
        Particle, PT, Simulation, Seeding
    ):
        label = create_label(particle, pT, simulation, seeding)
        run_single_particles(particle, pT, simulation, seeding, label)
