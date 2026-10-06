#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import uproot

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

    # Add the map writer(s)
    materialMapWriters = []
    # json map writer
    jsonOutputMapFormatsDict = {"json": JsonFormat.Json, "cbor": JsonFormat.Cbor}
    jsonOutputMapFormats = [
        jsonOutputMapFormatsDict[f]
        for f in outputMapFormats
        if f in jsonOutputMapFormatsDict
    ]
    if jsonOutputMapFormats:
        jmConverterCfg = MaterialMapJsonConverter.Config(
            processSensitives=True,
            processApproaches=True,
            processRepresenting=True,
            processBoundaries=True,
            processVolumes=False,
        )
        # Suffix for the map file is added in the writer depending on the format
        for writeFormat in jsonOutputMapFormats:
            materialMapWriters.append(
                JsonMaterialWriter(
                    level=loglevel,
                    converterCfg=jmConverterCfg,
                    fileName=outputFileBase + "_map",
                    writeFormat=writeFormat,
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


def plotEtaSurfaceDistance(
    filePath: str,
    treeName: str,
    outputSvg: str,
    etaRange: tuple[float, float] = (-4.0, 4.0),
    etaBins: int = 100,
    distRange: tuple[float, float] = (0.0, 1000.0),
    distBins: int = 100,
    excludeSurfaceIds: set[int] | None = None,
):
    # Distance between where material was actually found (Geant4 recording)
    # and the surface it got mapped to, as a function of the track eta.
    # Streamed and histogrammed in chunks: with O(1e6) tracks the full,
    # per-interaction (eta, distance) arrays are too large to hold in memory
    # at once.
    #
    # excludeSurfaceIds: raw GeometryIdentifier.value of surfaces to drop
    # (e.g. catch-all collector surfaces), matched exactly against sur_id.
    eta_edges = np.linspace(*etaRange, etaBins + 1)
    dist_edges = np.linspace(*distRange, distBins + 1)
    counts = np.zeros((etaBins, distBins))

    branches = ["v_eta", "sur_distance"]
    if excludeSurfaceIds:
        branches.append("sur_id")

    for batch in uproot.iterate(f"{filePath}:{treeName}", branches, library="np"):
        eta = batch["v_eta"]
        dist = batch["sur_distance"]  # jagged: one array per track
        nPerTrack = np.fromiter((len(d) for d in dist), dtype=np.int64, count=len(dist))
        eta_flat = np.repeat(eta, nPerTrack)
        dist_flat = np.concatenate(dist) if nPerTrack.sum() > 0 else np.empty(0)

        if excludeSurfaceIds:
            sur_id = batch["sur_id"]  # jagged, same shape as sur_distance
            id_flat = np.concatenate(sur_id) if nPerTrack.sum() > 0 else np.empty(0)
            keep = ~np.isin(id_flat, list(excludeSurfaceIds))
            eta_flat = eta_flat[keep]
            dist_flat = dist_flat[keep]

        h, _, _ = np.histogram2d(
            eta_flat, np.clip(dist_flat, *distRange), bins=[eta_edges, dist_edges]
        )
        counts += h

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(
        eta_edges,
        dist_edges,
        counts.T,
        norm=matplotlib.colors.LogNorm(vmin=1, vmax=max(counts.max(), 1)),
        cmap="viridis",
    )
    fig.colorbar(mesh, ax=ax, label="entries")
    ax.set_title("Material mapping distance vs. η")
    ax.set_xlabel("η")
    ax.set_ylabel("distance(true, mapped) [mm]")
    fig.savefig(outputSvg)
    plt.close(fig)


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
    p.add_argument(
        "--gen1",
        action="store_true",
        help="Map onto the Gen1 (Layer-based) geometry instead of Gen3 "
        "(default). Does not require re-recording: material recording is "
        "independent of the Gen1/Gen3 ACTS geometry split.",
    )

    args = p.parse_args()
    logLevel = logging.INFO
    gen3 = not args.gen1

    matDeco = None
    if args.matconfig != "":
        matDeco = acts.IMaterialDecorator.fromFile(args.matconfig)

    detector = getOpenDataDetector(matDeco, gen3=gen3)
    trackingGeometry = detector.trackingGeometry()

    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    runMaterialMapping(
        materialSurfaces,
        inputFile=Path(args.input),
        outputFileBase=args.output,
        outputMapFormats=["json", "root"],
        loglevel=logLevel,
        outputMaterialTracks=args.material_tracks_name,
        treeName=args.tree_name,
    ).run()

    # The calo material collector surfaces are just catch-alls to keep
    # calorimeter material from leaking onto tracker surfaces; they don't
    # represent a real, physically meaningful mapping location, so exclude
    # them from the distance diagnostic. Found by their geometric signature
    # (not TrackingGeometry.findVolumeByName: after portal fusion the
    # assigned GeometryIdentifier belongs to the fused portal, not to the
    # collector TrackingVolume's own geometryId). Constants mirror
    # addCaloMaterialCollector/addCaloMaterialCollectorDisc in
    # Plugins/DD4hep/src/OpenDataDetectorBuilder.cpp. Gen3-only: Gen1's
    # Layer-based geometry has no such collector volumes.
    excludeSurfaceIds = set()
    if gen3:
        for surface in materialSurfaces:
            bounds = surface.bounds
            if (
                isinstance(bounds, acts.CylinderBounds)
                and abs(bounds.values()[0] - 1215.0) < 0.5
                and abs(bounds.values()[1] - 3154.0) < 1.0
            ):
                excludeSurfaceIds.add(surface.geometryId.value)  # barrel collector
            elif (
                isinstance(bounds, acts.RadialBounds)
                and abs(bounds.values()[1] - 1225.0) < 0.5
            ):
                excludeSurfaceIds.add(surface.geometryId.value)  # disc collectors

    plotEtaSurfaceDistance(
        args.output + "_mapped.root",
        "material_tracks",  # RootMaterialTrackWriter's default tree name
        args.output + "_eta_distance.svg",
        excludeSurfaceIds=excludeSurfaceIds,
    )
