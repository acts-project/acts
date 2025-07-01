#!/usr/bin/env python3

import os
from pathlib import Path

import acts
from acts import MaterialMapJsonConverter
from acts.examples.odd import getOpenDataDetector
from acts.examples import (
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
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
            # if not os.path.isdir(outputDir + "/csv"):
            #    os.makedirs(outputDir + "/csv")
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(outputDir / "csv"),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO, outputDir=outputDir / "obj"
            )
            writer.write(context, trackingGeometry)

        if outputJson:
            # if not os.path.isdir(outputDir + "/json"):
            #    os.makedirs(outputDir + "/json")
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(outputDir / "json"),
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
                fileName=str(outputDir / "geometry-map"),
                writeFormat=JsonFormat.Json,
            )

            jmw.write(trackingGeometry)


if "__main__" == __name__:
    # detector = GenericDetector()
    detector = getOpenDataDetector()
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
