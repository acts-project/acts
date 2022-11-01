#!/usr/bin/env python3
from pathlib import Path

import acts
import acts.examples
from acts.examples.simulation import addParticleGun, addGeant4, EtaConfig
from acts.examples.odd import getOpenDataDetector
from acts.examples.geant4.dd4hep import DDG4DetectorConstruction
from common import getOpenDataDetectorDirectory

u = acts.UnitConstants


def runGeant4(
    g4detectorConstruction,
    trackingGeometry,
    field,
    outputDir,
    s: acts.examples.Sequencer = None,
):
    s = s or acts.examples.Sequencer(events=100, numThreads=1)
    s.config.logLevel = acts.logging.INFO
    rnd = acts.examples.RandomNumbers()
    addParticleGun(
        s,
        EtaConfig(-2.0, 2.0),
        rnd=rnd,
    )
    outputDir = Path(outputDir)
    addGeant4(
        s,
        g4detectorConstruction,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
    )
    return s


if "__main__" == __name__:
    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory()
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runGeant4(
        DDG4DetectorConstruction(detector), trackingGeometry, field, Path.cwd()
    ).run()
