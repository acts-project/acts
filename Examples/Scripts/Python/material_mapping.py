#!/usr/bin/env python3

import os

import acts
from acts import (
    SurfaceMaterialMapper,
    VolumeMaterialMapper,
    Navigator,
    Propagator,
    StraightLineStepper,
    MaterialMapJsonConverter,
)
from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
    MaterialMapping,
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
    mapSurface=True,
    mapVolume=True,
    readCachedSurfaceInformation=False,
    mappingStep=1,
    s=None,
):
    s = s or Sequencer(numThreads=1)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    wb = WhiteBoard(acts.logging.INFO)

    context = AlgorithmContext(0, 0, wb)

    for decorator in decorators:
        assert decorator.decorate(context) == ProcessCode.SUCCESS

    # Read material step information from a ROOT TTRee
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks="material-tracks",
            fileList=[
                os.path.join(
                    inputDir,
                    mapName + "_tracks.root"
                    if readCachedSurfaceInformation
                    else "geant4_material_tracks.root",
                )
            ],
            readCachedSurfaceInformation=readCachedSurfaceInformation,
        )
    )

    stepper = StraightLineStepper()

    mmAlgCfg = MaterialMapping.Config(context.geoContext, context.magFieldContext)
    mmAlgCfg.trackingGeometry = trackingGeometry
    mmAlgCfg.inputMaterialTracks = "material-tracks"

    if mapSurface:
        navigator = Navigator(
            trackingGeometry=trackingGeometry,
            resolveSensitive=True,
            resolveMaterial=True,
            resolvePassive=True,
        )
        propagator = Propagator(stepper, navigator)
        mapper = SurfaceMaterialMapper(level=acts.logging.INFO, propagator=propagator)
        mmAlgCfg.materialSurfaceMapper = mapper

    if mapVolume:
        navigator = Navigator(
            trackingGeometry=trackingGeometry,
        )
        propagator = Propagator(stepper, navigator)
        mapper = VolumeMaterialMapper(
            level=acts.logging.INFO, propagator=propagator, mappingStep=mappingStep
        )
        mmAlgCfg.materialVolumeMapper = mapper

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
        fileName=os.path.join(outputDir, mapName),
        writeFormat=JsonFormat.Json,
    )

    mmAlgCfg.materialWriters = [jmw]

    s.addAlgorithm(MaterialMapping(level=acts.logging.INFO, config=mmAlgCfg))

    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=mmAlgCfg.mappingMaterialCollection,
            filePath=os.path.join(
                outputDir,
                mapName + "_tracks.root",
            ),
            storeSurface=True,
            storeVolume=True,
        )
    )

    return s


if "__main__" == __name__:
    matDeco = acts.IMaterialDecorator.fromFile("geometry-map.json")
    detector, trackingGeometry, decorators = getOpenDataDetector(matDeco)

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=os.getcwd(),
        inputDir=os.getcwd(),
        readCachedSurfaceInformation=False,
    ).run()
