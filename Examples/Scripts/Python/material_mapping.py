#!/usr/bin/env python3

import argparse
from pathlib import Path
import acts

from acts import (
    Surface,
    MaterialMapper,
    IntersectionMaterialAssigner,
    BinnedSurfaceMaterialAccumulator,
    logging,
    GeometryContext,
)

from acts.json import MaterialMapJsonConverter

from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    MaterialMapping,
)

from acts.examples.root import (
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
    RootMaterialWriter,
)

from acts.examples.json import (
    JsonMaterialWriter,
    JsonFormat,
)

from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory


def runMaterialMapping(
    surfaces: list[Surface],
    inputFile: Path,
    outputFileBase: str,
    outputMapFormats: list[str] = ["json", "root"],
    loglevel: acts.logging.Level = acts.logging.INFO,
    outputMaterialTracks: str = "material_tracks",
    treeName: str = "material_tracks",
    readCachedSurfaceInformation: bool = False,
):
    # Create a sequencer
    print("Creating the sequencer with 1 thread (inter event information needed)")

    s = Sequencer(numThreads=1)

    # IO for material tracks reading
    wb = WhiteBoard(acts.logging.INFO)

    # Read material step information from a ROOT TTRee
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks=outputMaterialTracks,
            treeName=treeName,
            fileList=[str(inputFile)],
            readCachedSurfaceInformation=False,
        )
    )

    # Assignment setup : Intersection assigner
    materialAssingerConfig = IntersectionMaterialAssigner.Config()
    materialAssingerConfig.surfaces = surfaces
    materialAssinger = IntersectionMaterialAssigner(materialAssingerConfig, loglevel)

    # Accumulation setup : Binned surface material accumulater
    materialAccumulaterConfig = BinnedSurfaceMaterialAccumulator.Config()
    materialAccumulaterConfig.materialSurfaces = surfaces
    materialAccumulater = BinnedSurfaceMaterialAccumulator(
        materialAccumulaterConfig, loglevel
    )

    # Mapper setup
    materialMapperConfig = MaterialMapper.Config()
    materialMapperConfig.assignmentFinder = materialAssinger
    materialMapperConfig.surfaceMaterialAccumulater = materialAccumulater
    materialMapper = MaterialMapper(materialMapperConfig, loglevel)

    # Add the map writer(s)
    materialMapWriters = []
    # json map writer
    if "json" in outputMapFormats:
        jmConverterCfg = MaterialMapJsonConverter.Config(
            processSensitives=True,
            processApproaches=True,
            processRepresenting=True,
            processBoundaries=True,
            processVolumes=False,
        )
        # Suffix for the map file is added in the writer depending on the format
        materialMapWriters.append(
            JsonMaterialWriter(
                level=loglevel,
                converterCfg=jmConverterCfg,
                fileName=outputFileBase + "_map",
                writeFormat=JsonFormat.Json,
            )
        )
    if "root" in outputMapFormats:
        materialMapWriters.append(
            RootMaterialWriter(
                level=loglevel,
                filePath=outputFileBase + "_map.root",
            )
        )

    # Mapping Algorithm
    materialMappingConfig = MaterialMapping.Config()
    materialMappingConfig.materialMapper = materialMapper
    materialMappingConfig.inputMaterialTracks = outputMaterialTracks
    materialMappingConfig.mappedMaterialTracks = outputMaterialTracks + "_mapped"
    materialMappingConfig.unmappedMaterialTracks = outputMaterialTracks + "_unmapped"
    materialMappingConfig.materialWriters = materialMapWriters
    materialMapping = MaterialMapping(materialMappingConfig, loglevel)
    s.addAlgorithm(materialMapping)

    # Add the mapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=materialMappingConfig.mappedMaterialTracks,
            filePath=outputFileBase + "_mapped.root",
            storeSurface=True,
            storeVolume=False,
        )
    )

    # Add the unmapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=materialMappingConfig.unmappedMaterialTracks,
            filePath=outputFileBase + "_unmapped.root",
            storeSurface=True,
            storeVolume=False,
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
        "--matconfig", type=str, default="", help="Material configuration file"
    )

    args = p.parse_args()
    gContext = GeometryContext()
    logLevel = logging.INFO

    matDeco = None
    if args.matconfig != "":
        matDeco = acts.IMaterialDecorator.fromFile(args.matconfig)

    detector = getOpenDataDetector(matDeco)
    trackingGeometry = detector.trackingGeometry()

    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    runMaterialMapping(
        materialSurfaces,
        inputFile=Path(args.input),
        outputFileBase=args.output,
        outputMapFormats=["json", "root"],
        loglevel=logLevel,
    ).run()
