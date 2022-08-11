#!/usr/bin/env python3
import os
from acts.examples.odd import getOpenDataDetector
from acts.examples import (
    GenericDetector,
    AlignedDetector,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    SvgTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
)

import acts

from acts import MaterialMapJsonConverter


def runGeometry(
    trackingGeometry,
    decorators,
    outputDir,
    events=1,
    outputObj=True,
    outputCsv=True,
    outputJson=True,
    outputRoot=True,
    outputSvg=False,
):

    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0

        context = AlgorithmContext(ialg, ievt, eventStore)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            csvOutputDir = os.path.join(outputDir, "csv")
            if not os.path.exists(csvOutputDir):
                os.makedirs(csvOutputDir)
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=csvOutputDir,
                writePerEvent=True,
            )

            writer.write(context)

        if outputObj:
            outputDirObj = os.path.join(outputDir, "obj")
            if not os.path.exists(outputDirObj):
                os.makedirs(outputDirObj)
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO, outputDir=outputDirObj
            )
            writer.write(context, trackingGeometry)

        if outputSvg:
            outputDirSvg = os.path.join(outputDir, "svg")
            if not os.path.exists(outputDirSvg):
                os.makedirs(outputDirSvg)
            svgWriterCfg = SvgTrackingGeometryWriter.Config(outputDir=outputDirSvg)
            svgWriter = SvgTrackingGeometryWriter(
                level=acts.logging.INFO, config=svgWriterCfg
            )
            svgWriter.write(context, trackingGeometry)

        if outputJson:
            outputDirJson = os.path.join(outputDir, "json")
            if not os.path.exists(outputDirJson):
                os.makedirs(outputDirJson)
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=outputDirJson,
                writePerEvent=True,
                writeSensitive=True,
            )
            writer.write(context)

            jmConverterCfg = MaterialMapJsonConverter.Config(
                processSensitives=True,
                processApproaches=True,
                processRepresenting=True,
                processBoundaries=True,
                processVolumes=True,
                context=context.geoContext,
            )

            jmw = JsonMaterialWriter(
                level=acts.logging.VERBOSE,
                converterCfg=jmConverterCfg,
                fileName=os.path.join(outputDirJson, "geometry-map"),
                writeFormat=JsonFormat.Json,
            )

            jmw.write(trackingGeometry)


if "__main__" == __name__:
    # detector, trackingGeometry, decorators = AlignedDetector.create()
    detector, trackingGeometry, decorators = GenericDetector.create()
    # detector, trackingGeometry, decorators = getOpenDataDetector(getOpenDataDetectorDirectory() )

    runGeometry(trackingGeometry, decorators, outputDir=os.getcwd())
