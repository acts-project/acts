#!/bin/python3
import sys
from pathlib import Path
import os

import acts
import acts.examples

from truth_tracking_gsf import runTruthTrackingGsf

u = acts.UnitConstants

if __name__ == "__main__":
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    digiConfigFile = (
        srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    )
    assert digiConfigFile.exists()

    detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    seq = acts.examples.Sequencer(
        events=1,
        numThreads=1,
        trackFpes=False,
        logLevel=acts.logging.VERBOSE,
    )

    seq.addWriter(
        acts.examples.CsvTrackingGeometryWriter(
            level=acts.logging.INFO,
            trackingGeometry=trackingGeometry,
            outputDir=Path.cwd(),
            writePerEvent=False,
        )
    )

    runTruthTrackingGsf(
        s=seq,
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        outputDir=Path.cwd(),
    ).run()
