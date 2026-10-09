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
)

from acts.examples import (
    Sequencer,
    MaterialMapping,
)

from acts.examples.root import (
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
    RootMaterialWriter,
)

from acts.json import TrackingGeometryMaterialJsonConverter

from acts.examples.odd import getOpenDataDetector


def runMaterialMapping(
    surfaces: list[Surface],
    inputFile: Path,
    outputFileBase: str,
    loglevel: acts.logging.Level = acts.logging.INFO,
    outputMaterialTracks: str = "material_tracks",
    treeName: str = "material_tracks",
    s: Sequencer | None = None,
):
    """Configure mapping and return the sequencer and algorithm.

    After running the sequencer, serialize the algorithm's ``material`` result.
    """
    if s is None:
        s = Sequencer(numThreads=1)

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

    # Accumulation setup : Binned surface material accumulator
    materialAccumulatorConfig = BinnedSurfaceMaterialAccumulator.Config()
    materialAccumulatorConfig.materialSurfaces = surfaces
    materialAccumulator = BinnedSurfaceMaterialAccumulator(
        materialAccumulatorConfig, loglevel
    )

    # Mapper setup
    materialMapperConfig = MaterialMapper.Config()
    materialMapperConfig.assignmentFinder = materialAssinger
    materialMapperConfig.surfaceMaterialAccumulator = materialAccumulator
    materialMapper = MaterialMapper(materialMapperConfig, loglevel)

    # Mapping Algorithm
    materialMappingConfig = MaterialMapping.Config()
    materialMappingConfig.materialMapper = materialMapper
    materialMappingConfig.inputMaterialTracks = outputMaterialTracks
    materialMappingConfig.mappedMaterialTracks = outputMaterialTracks + "_mapped"
    materialMappingConfig.unmappedMaterialTracks = outputMaterialTracks + "_unmapped"
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

    return s, materialMapping


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

    p.add_argument(
        "--tree-name",
        type=str,
        default="material_tracks",
        help="Input material track tree name",
    )

    p.add_argument(
        "--material_tracks-name",
        type=str,
        default="material_tracks",
        help="Input material track collection name",
    )

    args = p.parse_args()
    logLevel = logging.INFO

    matDeco = None
    if args.matconfig != "":
        matDeco = acts.IMaterialDecorator.fromFile(args.matconfig)

    detector = getOpenDataDetector(matDeco)
    trackingGeometry = detector.trackingGeometry()

    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    s, mapping = runMaterialMapping(
        materialSurfaces,
        inputFile=Path(args.input),
        outputFileBase=args.output,
        loglevel=logLevel,
        outputMaterialTracks=args.material_tracks_name,
        treeName=args.tree_name,
        s=Sequencer(events=args.events, numThreads=1),
    )
    s.run()
    TrackingGeometryMaterialJsonConverter().toFile(
        mapping.material, args.output + "_map.json"
    )
    RootMaterialWriter(
        level=logLevel, filePath=args.output + "_map.root"
    ).writeMaterial(mapping.material)
