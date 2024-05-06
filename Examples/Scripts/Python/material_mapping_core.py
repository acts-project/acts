#!/usr/bin/env python3

import argparse

import acts
from acts import (
    MaterialMapper,
    IntersectionMaterialAssigner,
    BinnedSurfaceMaterialAccumulater,
    MaterialMapJsonConverter,
)

from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
    CoreMaterialMapping,
    JsonMaterialWriter,
    RootMaterialWriter,
    JsonFormat,
)

from acts.examples.dd4hep import (
    DD4hepDetector,
    DD4hepDetectorOptions,
    DD4hepGeometryService,
)

from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory


def runMaterialMapping(surfaces, inputFile, outputFile, outputMap, loglevel):
    # Create a sequencer
    print("Creating the sequencer with 1 thread (inter event information needed)")

    s = Sequencer(numThreads=1)

    # IO for material tracks reading
    wb = WhiteBoard(acts.logging.INFO)

    # Read material step information from a ROOT TTRee
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks="material-tracks",
            fileList=[inputFile],
            readCachedSurfaceInformation=False,
        )
    )

    # Assignment setup : Intersection assigner
    materialAssingerConfig = IntersectionMaterialAssigner.Config()
    materialAssingerConfig.surfaces = surfaces
    materialAssinger = IntersectionMaterialAssigner(materialAssingerConfig, loglevel)

    # Accumulation setup : Binned surface material accumulater
    materialAccumulaterConfig = BinnedSurfaceMaterialAccumulater.Config()
    materialAccumulaterConfig.materialSurfaces = surfaces
    materialAccumulater = BinnedSurfaceMaterialAccumulater(
        materialAccumulaterConfig, loglevel
    )

    # Mapper setup
    materialMapperConfig = MaterialMapper.Config()
    materialMapperConfig.assignmentFinder = materialAssinger
    materialMapperConfig.surfaceMaterialAccumulater = materialAccumulater
    materialMapper = MaterialMapper(materialMapperConfig, loglevel)

    # Add the map writer(s)
    mapWriters = []
    # json map writer
    context = AlgorithmContext(0, 0, wb)
    jmConverterCfg = MaterialMapJsonConverter.Config(
        processSensitives=True,
        processApproaches=True,
        processRepresenting=True,
        processBoundaries=True,
        processVolumes=True,
        context=context.geoContext,
    )
    mapWriters.append(
        JsonMaterialWriter(
            level=loglevel,
            converterCfg=jmConverterCfg,
            fileName=outputMap + "",
            writeFormat=JsonFormat.Json,
        )
    )
    mapWriters.append(RootMaterialWriter(level=loglevel, filePath=outputMap + ".root"))

    # Mapping Algorithm
    coreMaterialMappingConfig = CoreMaterialMapping.Config()
    coreMaterialMappingConfig.materialMapper = materialMapper
    coreMaterialMappingConfig.inputMaterialTracks = "material-tracks"
    coreMaterialMappingConfig.mappedMaterialTracks = "mapped-material-tracks"
    coreMaterialMappingConfig.unmappedMaterialTracks = "unmapped-material-tracks"
    coreMaterialMappingConfig.materiaMaplWriters = mapWriters
    coreMaterialMapping = CoreMaterialMapping(coreMaterialMappingConfig, loglevel)
    s.addAlgorithm(coreMaterialMapping)

    # Add the mapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks="mapped-material-tracks",
            filePath=outputFile + "_mapped.root",
            storeSurface=True,
            storeVolume=True,
        )
    )

    # Add the unmapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks="unmapped-material-tracks",
            filePath=outputFile + "_unmapped.root",
            storeSurface=True,
            storeVolume=True,
        )
    )

    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser()

    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to process"
    )
    p.add_argument(
        "-i", "--input", type=str, default="", help="Input file with material tracks"
    )
    p.add_argument(
        "-o", "--output", type=str, default="", help="Output file (core) name"
    )
    p.add_argument(
        "-m", "--map", type=str, default="", help="Output file for the material map"
    )

    p.add_argument(
        "--matconfig", type=str, default="", help="Material configuration file"
    )

    p.add_argument(
        "--experimental",
        action=argparse.BooleanOptionalAction,
        help="Construct experimental geometry",
    )

    args = p.parse_args()

    if args.experimental:
        odd_xml = getOpenDataDetectorDirectory() / "xml" / "OpenDataDetector.xml"

        # Create the dd4hep geometry service and detector
        dd4hepConfig = DD4hepGeometryService.Config()
        dd4hepConfig.logLevel = acts.logging.INFO
        dd4hepConfig.xmlFileNames = [str(odd_xml)]
        dd4hepGeometryService = DD4hepGeometryService(dd4hepConfig)
        dd4hepDetector = DD4hepDetector(dd4hepGeometryService)

        cOptions = DD4hepDetectorOptions(logLevel=acts.logging.INFO, emulateToGraph="")

        # Context and options
        geoContext = acts.GeometryContext()
        [detector, contextors, store] = dd4hepDetector.finalize(geoContext, cOptions)

        materialSurfaces = detector.extractMaterialSurfaces()

    else:
        matDeco = None
        if args.matconfig != "":
            matDeco = acts.IMaterialDecorator.fromFile(args.matconfig)

        [detector, trackingGeometry, decorators] = getOpenDataDetector(matDeco)

        materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    runMaterialMapping(
        materialSurfaces, args.input, args.output, args.map, acts.logging.INFO
    ).run()
