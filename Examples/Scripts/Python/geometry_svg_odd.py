#!/usr/bin/env python3
import os
from acts.examples.odd import getOpenDataDetector
from acts.examples import (
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    GenericDetector,
    SvgTrackingGeometryWriter
)
from pathlib import Path
import acts


def runGeometry(
        trackingGeometry,
        outputDir,
        events=1):

    for ievt in range(events):
        eventStore = WhiteBoard(
            name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0

        context = AlgorithmContext(ialg, ievt, eventStore)

        outputDirSvg = os.path.join(outputDir, "svg")
        if not os.path.exists(outputDirSvg):
            os.makedirs(outputDirSvg)
        svgWriterCfg = SvgTrackingGeometryWriter.Config(outputDir=outputDirSvg)
        svgWriter = SvgTrackingGeometryWriter(
            level=acts.logging.INFO, config=svgWriterCfg)
        svgWriter.write(context, trackingGeometry)


if "__main__" == __name__:
    detector, trackingGeometry, decorators = GenericDetector.create()
    
    runGeometry(trackingGeometry, outputDir=os.getcwd())
