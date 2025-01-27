#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples
from acts.examples.simulation import addParticleGun, addGeant4, EtaConfig
from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants


def runGeant4(
    detector,
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
        detector,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
    )
    return s


if "__main__" == __name__:
    detector, trackingGeometry, decorators = getOpenDataDetector()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runGeant4(detector, trackingGeometry, field, Path.cwd()).run()
