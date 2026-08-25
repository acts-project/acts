"""ODD KF alignment workflow — simulation stage: particle gun, Fatras, digitization."""

from __future__ import annotations

from pathlib import Path

import acts
import acts.examples

from kfalignment_config import DEFAULT_NUM_EVENTS, DEFAULT_NUM_THREADS

u = acts.UnitConstants


def runSimulation(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    detector=None,
    numEvents: int = DEFAULT_NUM_EVENTS,
    s: acts.examples.Sequencer = None,
):
    """Run particle gun + Fatras + digitization and write CSV/ROOT outputs."""
    from acts.examples.simulation import (
        addParticleGun,
        ParticleConfig,
        EtaConfig,
        PhiConfig,
        MomentumConfig,
        addFatras,
        addDigitization,
    )

    s = s or acts.examples.Sequencer(
        events=numEvents, numThreads=DEFAULT_NUM_THREADS, logLevel=acts.logging.INFO
    )

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)
    outputDir.mkdir(exist_ok=True)

    logger = acts.getDefaultLogger("Simulation", acts.logging.INFO)

    addParticleGun(
        s,
        ParticleConfig(num=1, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
        EtaConfig(-2.0, 2.0, uniform=True),
        MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
        PhiConfig(0.0 * u.degree, 360.0 * u.degree),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(2.0 * u.mm, 2.0 * u.mm, 2.0 * u.mm, 2.0 * u.ns),
        ),
        multiplicity=1,
        rnd=rnd,
        outputDirCsv=outputDir,
        outputDirRoot=outputDir,
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        enableInteractions=True,
        outputDirCsv=outputDir,
        outputDirRoot=outputDir,
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
        outputDirCsv=outputDir,
        outputDirRoot=outputDir,
    )

    logger.info("Simulation output will be saved to: {}", str(outputDir))

    return s
