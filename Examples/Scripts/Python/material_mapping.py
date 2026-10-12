#!/usr/bin/env python3

import argparse
from pathlib import Path
import acts

from acts import (
    Surface,
    MaterialMapper,
    IntersectionMaterialAssigner,
    BinnedSurfaceMaterialAccumulator,
    GridSurfaceMaterialAccumulator,
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

from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory


def runMaterialMapping(
    surfaces: list[Surface],
    inputFile: Path,
    outputFileBase: str,
    loglevel: acts.logging.Level = acts.logging.INFO,
    outputMaterialTracks: str = "material_tracks",
    treeName: str = "material_tracks",
    s: Sequencer | None = None,
    accumulator: str = "binned",
):
    """Configure binned or grid mapping and return the sequencer and algorithm.

    After running the sequencer, serialize the algorithm's ``material`` result.
    """
    accumulators = {
        "binned": BinnedSurfaceMaterialAccumulator,
        "grid": GridSurfaceMaterialAccumulator,
    }
    if accumulator not in accumulators:
        raise ValueError("accumulator must be 'binned' or 'grid'")
    accumulatorClass = accumulators[accumulator]

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

    # Accumulation setup
    materialAccumulatorConfig = accumulatorClass.Config()
    materialAccumulatorConfig.materialSurfaces = surfaces
    materialAccumulator = accumulatorClass(materialAccumulatorConfig, loglevel)

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


def main():
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
        "--matconfig",
        type=str,
        default="",
        help="Material configuration file (grid default: ODD mapping configuration)",
    )

    p.add_argument(
        "--accumulator",
        choices=["binned", "grid"],
        default="binned",
        help="Material accumulator to use (default: binned)",
    )
    p.add_argument(
        "--output-map-formats",
        choices=["json", "cbor", "root"],
        nargs="+",
        help="Map formats (default: json for grid; json and root for binned)",
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
    outputMapFormats = args.output_map_formats
    if outputMapFormats is None:
        outputMapFormats = ["json"] if args.accumulator == "grid" else ["json", "root"]
    if args.accumulator == "grid" and "root" in outputMapFormats:
        p.error("Grid material maps require JSON or CBOR output")

    matconfig = args.matconfig
    if matconfig == "" and args.accumulator == "grid":
        matconfig = (
            getOpenDataDetectorDirectory()
            / "config"
            / "odd-material-mapping-config.json"
        )
        print(f"Using ODD prototype material configuration: {matconfig}")

    matDeco = None
    if matconfig != "":
        matDeco = acts.IMaterialDecorator.fromFile(matconfig)

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
        accumulator=args.accumulator,
    )
    s.run()
    for outputFormat in outputMapFormats:
        if outputFormat == "root":
            RootMaterialWriter(
                level=logLevel, filePath=args.output + "_map.root"
            ).writeMaterial(mapping.material)
        else:
            TrackingGeometryMaterialJsonConverter().toFile(
                mapping.material, args.output + "_map." + outputFormat
            )


if "__main__" == __name__:
    main()
