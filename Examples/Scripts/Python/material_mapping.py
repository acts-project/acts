#!/usr/bin/env python3

import os
import argparse

import acts
from acts import (
    MaterialMapper,
    IntersectionMaterialAssigner,
    BinnedSurfaceMaterialAccumulater,
)

from acts.json import MaterialMapJsonConverter

from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    MaterialMapping,
)

from acts.examples.root import (
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
)

from acts.examples.json import (
    JsonMaterialWriter,
    JsonFormat,
)

from acts.examples.odd import getOpenDataDetector


def runMaterialMapping(
    trackingGeometry,
    decorators,
    outputDir,
    inputDir,
    mapName="material-map",
    mapFormat=JsonFormat.Json,
    mapSurface=True,
    mapVolume=False,
    readCachedSurfaceInformation=False,
    mappingStep=1,
    s=None,
):
    del mapVolume
    del mappingStep

    if not mapSurface:
        raise ValueError("Surface mapping must be enabled.")

    s = s or Sequencer(numThreads=1)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    wb = WhiteBoard(acts.logging.INFO)
    context = AlgorithmContext(0, 0, wb, 0)

    for decorator in decorators:
        assert decorator.decorate(context) == ProcessCode.SUCCESS

    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks="material-tracks",
            fileList=[
                os.path.join(
                    inputDir,
                    (
                        mapName + "_tracks.root"
                        if readCachedSurfaceInformation
                        else "geant4_material_tracks.root"
                    ),
                )
            ],
            readCachedSurfaceInformation=readCachedSurfaceInformation,
        )
    )

    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    assignmentCfg = IntersectionMaterialAssigner.Config()
    assignmentCfg.surfaces = materialSurfaces
    assignmentFinder = IntersectionMaterialAssigner(assignmentCfg, acts.logging.INFO)

    accumulaterCfg = BinnedSurfaceMaterialAccumulater.Config()
    accumulaterCfg.materialSurfaces = materialSurfaces
    surfaceAccumulater = BinnedSurfaceMaterialAccumulater(
        accumulaterCfg, acts.logging.INFO
    )

    mapperCfg = MaterialMapper.Config()
    mapperCfg.assignmentFinder = assignmentFinder
    mapperCfg.surfaceMaterialAccumulater = surfaceAccumulater
    materialMapper = MaterialMapper(mapperCfg, acts.logging.INFO)

    mmAlgCfg = MaterialMapping.Config(context.geoContext)
    mmAlgCfg.inputMaterialTracks = "material-tracks"
    mmAlgCfg.mappedMaterialTracks = "mapped-material-tracks"
    mmAlgCfg.unmappedMaterialTracks = "unmapped-material-tracks"
    mmAlgCfg.materialMapper = materialMapper

    jmConverterCfg = MaterialMapJsonConverter.Config(
        processSensitives=True,
        processApproaches=True,
        processRepresenting=True,
        processBoundaries=True,
        processVolumes=False,
        context=context.geoContext,
    )

    jmw = JsonMaterialWriter(
        level=acts.logging.VERBOSE,
        converterCfg=jmConverterCfg,
        fileName=os.path.join(outputDir, mapName),
        writeFormat=mapFormat,
    )

    mmAlgCfg.materialWriters = [jmw]

    s.addAlgorithm(MaterialMapping(level=acts.logging.INFO, config=mmAlgCfg))

    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=mmAlgCfg.mappedMaterialTracks,
            filePath=os.path.join(
                outputDir,
                mapName + "_tracks.root",
            ),
            storeSurface=True,
            storeVolume=False,
        )
    )

    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser(description="Script to generate ACTS material map")
    p.add_argument(
        "-o",
        "--outFile",
        type=str,
        default="material-map.json",
        help="Output filename for the generated material map. Supported formats: JSON, CBOR.",
    )
    args = p.parse_args()
    if ".json" in args.outFile:
        mapFormat = JsonFormat.Json
    elif ".cbor" in args.outFile:
        mapFormat = JsonFormat.Cbor
    else:
        print(
            "ERROR(material_mapping.py): please provide an output name ending with .json or .cbor"
        )
        exit()

    mapName = args.outFile.split(".")[0]

    detector = getOpenDataDetector(None)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=os.getcwd(),
        inputDir=os.getcwd(),
        readCachedSurfaceInformation=False,
        mapName=mapName,
        mapFormat=mapFormat,
    ).run()
