#!/usr/bin/env python3

import os
import json
from pathlib import Path

import acts
from acts import MaterialMapJsonConverter
from acts.examples.odd import getOpenDataDetector
import argparse
from acts.examples import (
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
    GenericDetector,
)


def runGeometry(
    trackingGeometry,
    decorators,
    outputDir: Path,
    events=1,
    outputObj=True,
    outputCsv=True,
    outputJson=True,
):
    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0
        ithread = 0

        context = AlgorithmContext(ialg, ievt, eventStore, ithread)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            csvDir = outputDir / "csv"
            if not csvDir.exists():
                csvDir.mkdir(parents=True)
            elif not csvDir.is_dir():
                raise RuntimeError(f"Output directory {csvDir} is not a directory")
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(csvDir),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            objDir = outputDir / "obj"
            if not objDir.exists():
                objDir.mkdir(parents=True)
            elif not objDir.is_dir():
                raise RuntimeError(f"Output directory {objDir} is not a directory")
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO, outputDir=str(objDir)
            )
            writer.write(context, trackingGeometry)

        if outputJson:
            jsonDir = outputDir / "json"
            if not jsonDir.exists():
                jsonDir.mkdir(parents=True)
            elif not jsonDir.is_dir():
                raise RuntimeError(f"Output directory {jsonDir} is not a directory")
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(jsonDir),
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
                processNonMaterial=True,
                context=context.geoContext,
            )

            jmw = JsonMaterialWriter(
                level=acts.logging.VERBOSE,
                converterCfg=jmConverterCfg,
                fileName=str(jsonDir / "geometry-map"),
                writeFormat=JsonFormat.Json,
            )

            jmw.write(trackingGeometry)


if "__main__" == __name__:
    p = argparse.ArgumentParser()
    p.add_argument("--verbosity", "-v", default=0, action="count")
    args = p.parse_args()
    logLevel = acts.logging.INFO
    if args.verbosity == 1:
        logLevel = acts.logging.DEBUG
    elif args.verbosity >= 2:
        logLevel = acts.logging.VERBOSE

    # detector = AlignedDetector()
    detector = GenericDetector(gen3=True, logLevel=logLevel)
    # detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    runGeometry(trackingGeometry, decorators, outputDir=Path.cwd())

    # Uncomment if you want to create the geometry id mapping for DD4hep
    # dd4hepIdGeoIdMap = acts.examples.dd4hep.createDD4hepIdGeoIdMap(trackingGeometry)
    # dd4hepIdGeoIdValueMap = {}
    # for key, value in dd4hepIdGeoIdMap.items():
    #     dd4hepIdGeoIdValueMap[key] = value.value()

    # with open('odd-dd4hep-geoid-mapping.json', 'w') as outfile:
    #    json.dump(dd4hepIdGeoIdValueMap, outfile)
