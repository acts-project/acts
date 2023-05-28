#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil
from enum import Enum
import itertools
import contextlib

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
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    addCKFTracks,
)

from physmon_common import makeSetup

u = acts.UnitConstants

setup = makeSetup()

Simulation = Enum("Simulation", ["Fatras", "Geant4"])


def run_single_particles(particle, pT, simulation, label):
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(events=10, numThreads=1)

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
                setup.detector,
                setup.trackingGeometry,
                setup.field,
                rnd=rnd,
            )
        else:
            raise ValueError(f"unhandled simulation: {simulation}")

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
            SeedFinderOptionsArg(bFieldInZ=2 * u.T, beamPos=(0.0, 0.0)),
            seedingAlgorithm=SeedingAlgorithm.Default,
            geoSelectionConfigFile=setup.geoSel,
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


def create_label(particle, pT, simulation):
    return (
        f"{acts.pdgToShortAbsString(particle)}_{int(pT)}GeV_{simulation.name.lower()}"
    )


# disabled FPE monitoring for now because of G4
# TODO use acts.FpeMonitor()
with contextlib.nullcontext():
    for particle, pt, simulation in itertools.product(
        [
            acts.PdgParticle.eMuon,
            acts.PdgParticle.ePionPlus,
            acts.PdgParticle.eElectron,
        ],
        [1 * u.GeV, 10 * u.GeV, 100 * u.GeV],
        [Simulation.Fatras],  # TODO Simulation.Geant4
    ):
        label = create_label(particle, pt, simulation)
        run_single_particles(particle, pt, simulation, label)
