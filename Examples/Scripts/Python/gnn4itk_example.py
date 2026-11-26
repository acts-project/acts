#!/usr/bin/env python3

"""
Example: GNN track finding with module maps on ATLAS ITk data.

Demonstrates reading pre-simulated ATLAS data and running GNN with module map
(geometry-based) graph construction. Typical workflow for ATLAS ITk.

All parameters are hardcoded except file paths. Users should copy and modify this script
for their specific needs.
"""

from pathlib import Path
import argparse

import acts
import acts.examples
import acts.examples.root
from acts.examples.reconstruction import (
    addGnn,
    addSpacePointsMaking,
    addTrackWriters,
)
from acts.examples.simulation import ParticleSelectorConfig, addDigiParticleSelection
from acts.gnn import (
    ModuleMapCuda,
    CudaTrackBuilding,
    DWalkTrackBuilding,
    EdgeLayerConnector,
)
from acts.examples.gnn import NodeFeature

u = acts.UnitConstants


def runGNN4ITk(
    inputRootDump: Path,
    moduleMapPath: str,
    gnnModel: Path,
    outputDir: Path = Path.cwd(),
    events: int = 1,
    bufferEvents: int | None = None,
    useEdgeLayerConnector: bool = False,
    trackBuilder: str = "junctionremoval",
    edgeClassifierCut: float | None = None,
    pixelOnly: bool = False,
    numFeatures: int = 12,
    trackingGeometry=None,
    field=None,
    geometryIdMap=None,
    geoSelectionConfigFile: Path | None = None,
    logLevel=acts.logging.INFO,
):
    """
    Run GNN tracking with module maps on ATLAS Athena dumps.

    This example shows reading pre-simulated ATLAS data and running GNN with
    geometry-based graph construction. All GNN parameters hardcoded.

    Args:
        inputRootDump: Path to input ROOT file (ATLAS Athena dump format)
        moduleMapPath: Path prefix for module map files
                      (will load .doublets.root and .triplets.root)
        gnnModel: Path to trained model (.pt, .onnx, or .engine)
        outputDir: Output directory for performance files
        events: Number of events to process
        useEdgeLayerConnector: Use EdgeLayerConnector instead of CudaTrackBuilding
        trackBuilder: Track builder to run: junctionremoval, edgelayerconnector, or dwalk
        edgeClassifierCut: Optional edge classifier cut. Defaults to 0.01 for D-WALK and
                           0.5 for the other builders.
        pixelOnly: Run on pixel space points only (for pixel-only GNN models)
        numFeatures: Number of node features the model expects: 4 (space point
                     R, Phi, Z, Eta) or 12 (space point + 2 clusters)
        trackingGeometry: Optional ITk tracking geometry. If given, the measurements
                          are projected onto its surfaces and the CKF is run from
                          the GNN tracks, e.g. to extend pixel-only tracks into
                          the strips
        field: Magnetic field, needed for the CKF
        geometryIdMap: Athena <-> ACTS geometry identifier map, needed for the CKF
                       (see acts.examples.readCsvGeometryIdMap)
        geoSelectionConfigFile: JSON file with the pixel volumes of the tracking
                                geometry, needed for the CKF seeding space points
        logLevel: Logging level
    """
    # Validate inputs
    assert inputRootDump.exists(), f"Input file not found: {inputRootDump}"
    assert Path(
        moduleMapPath + ".doublets.root"
    ).exists(), f"Module map not found: {moduleMapPath}.doublets.root"
    assert Path(
        moduleMapPath + ".triplets.root"
    ).exists(), f"Module map not found: {moduleMapPath}.triplets.root"
    assert gnnModel.exists(), f"Model file not found: {gnnModel}"
    if numFeatures not in (4, 12):
        raise ValueError(f"Unsupported number of features: {numFeatures}")
    if trackingGeometry is not None and (
        field is None or geometryIdMap is None or geoSelectionConfigFile is None
    ):
        raise ValueError(
            "The CKF needs a magnetic field, a geometry id map and a geometry "
            "selection file"
        )

    # Avoid an expensive recursive source-tree scan for automatic FPE masks on
    # large CFS-mounted ACTS workspaces.
    acts.examples.Sequencer._autoFpeMasks = []
    s = acts.examples.Sequencer(
        events=events,
        numThreads=1,
    )

    # Read ATLAS Athena ROOT dump
    reader = acts.examples.root.RootAthenaDumpReader(
        level=logLevel,
        treename="GNN4ITk",
        inputfiles=[str(inputRootDump)],
        outputSpacePoints="spacepoints",
        outputPixelSpacePoints="pixel_spacepoints",
        outputClusters="clusters",
        outputMeasurements="measurements",
        outputMeasurementSubset="measurement_subset",
        outputMeasurementParticlesMap="measurement_particles_map",
        outputParticleMeasurementsMap="particle_measurements_map",
        outputParticles="particles",
        skipOverlapSPsPhi=True,
        skipOverlapSPsEta=False,
        absBoundaryTolerance=0.01 * u.mm,
        trackingGeometry=trackingGeometry,
        geometryIdMap=geometryIdMap,
    )

    if bufferEvents is not None:
        s.addReader(
            acts.examples.BufferedReader(
                level=logLevel,
                upstreamReader=reader,
                bufferSize=min(bufferEvents, events),
            )
        )
    else:
        s.addReader(reader)

    # Select primary particles with minimum 7 hits and 1 GeV pT for efficiency evaluation
    s.addWhiteboardAlias("particles_simulated_selected", "particles")
    particleSelectorConfig = ParticleSelectorConfig(
        pt=(1.0 * u.GeV, None),
        hits=(7, None),
        removeSecondaries=True,
    )
    addDigiParticleSelection(s, particleSelectorConfig, logLevel=logLevel)

    # Configure GNN stages for module map workflow
    # All parameters hardcoded based on ITk configuration

    # Stage 1: Graph construction via module map
    moduleMapConfig = {
        "level": logLevel,
        "moduleMapPath": moduleMapPath,
        "rScale": 1000.0,
        "phiScale": 3.141592654,
        "zScale": 1000.0,
        "etaScale": 1.0,
        "gpuDevice": 0,
        "gpuBlocks": 512,
    }
    graphConstructor = ModuleMapCuda(**moduleMapConfig)

    if useEdgeLayerConnector:
        trackBuilder = "edgelayerconnector"
    trackBuilder = trackBuilder.lower()
    if trackBuilder not in {"junctionremoval", "edgelayerconnector", "dwalk"}:
        raise ValueError(f"Unsupported track builder: {trackBuilder}")

    defaultEdgeClassifierCut = 0.01 if trackBuilder == "dwalk" else 0.5
    resolvedEdgeClassifierCut = (
        defaultEdgeClassifierCut if edgeClassifierCut is None else edgeClassifierCut
    )

    # Stage 2: Single-stage edge classification (auto-detect backend)
    gnnModel = Path(gnnModel)
    edgeClassifierConfig = {
        "level": logLevel,
        "modelPath": str(gnnModel),
        "cut": resolvedEdgeClassifierCut,
    }

    if gnnModel.suffix == ".pt":
        edgeClassifierConfig["useEdgeFeatures"] = True
        from acts.gnn import TorchEdgeClassifier

        edgeClassifiers = [TorchEdgeClassifier(**edgeClassifierConfig)]
    elif gnnModel.suffix == ".onnx":
        from acts.gnn import OnnxEdgeClassifier

        edgeClassifiers = [OnnxEdgeClassifier(**edgeClassifierConfig)]
    elif gnnModel.suffix == ".engine":
        from acts.gnn import TensorRTEdgeClassifier

        edgeClassifiers = [TensorRTEdgeClassifier(**edgeClassifierConfig)]
    else:
        raise ValueError(f"Unsupported model format: {gnnModel.suffix}")

    # Stage 3: Track building
    if trackBuilder == "edgelayerconnector":
        trackBuilderObj = EdgeLayerConnector(
            level=logLevel,
            blockSize=512,
            maxHitsPerTrack=30,
            minHits=3,
            weightsCut=0.01,
        )
    elif trackBuilder == "dwalk":
        trackBuilderObj = DWalkTrackBuilding(
            level=logLevel,
            thMin=0.1,
            thAdd=0.6,
            minCandidateSize=3,
            radialFeatureIndex=0,
            pathMetric="score_weighted_length",
        )
    else:
        trackBuilderObj = CudaTrackBuilding(
            level=logLevel,
            useOneBlockImplementation=False,
            doJunctionRemoval=True,
        )

    # Node features: ITk 12-feature configuration (space point + 2 clusters),
    # truncated to the space point features for 4-feature models
    e = NodeFeature
    nodeFeatures = [
        e.R,
        e.Phi,
        e.Z,
        e.Eta,
        e.Cluster1R,
        e.Cluster1Phi,
        e.Cluster1Z,
        e.Cluster1Eta,
        e.Cluster2R,
        e.Cluster2Phi,
        e.Cluster2Z,
        e.Cluster2Eta,
    ]
    featureScales = [1000.0, 3.141592654, 1000.0, 1.0] * 3
    nodeFeatures = nodeFeatures[:numFeatures]
    featureScales = featureScales[:numFeatures]

    spacePoints = "pixel_spacepoints" if pixelOnly else "spacepoints"

    # Add GNN tracking (space points already created by reader, so no trackingGeometry)
    addGnn(
        s,
        graphConstructor=graphConstructor,
        edgeClassifiers=edgeClassifiers,
        trackBuilder=trackBuilderObj,
        nodeFeatures=nodeFeatures,
        featureScales=featureScales,
        inputSpacePoints=spacePoints,
        inputClusters="clusters",
        outputDirRoot=str(outputDir),
        geometryIdMap=geometryIdMap,
        logLevel=logLevel,
    )

    # Stage 4 (optional): CKF from the GNN tracks, e.g. into the strips
    if trackingGeometry is not None:
        # Space points for the seeds made from the measurements, as in the
        # standard seeding
        s.addWhiteboardAlias("ckf_measurement_subset", "measurement_subset")
        ckfSpacePoints = addSpacePointsMaking(
            s,
            trackingGeometry,
            geoSelectionConfigFile,
            None,
            logLevel,
            prefix="ckf_",
        )

        parametersConfig = {
            "level": logLevel,
            "inputProtoTracks": "gnn-protoTracks",
            "inputSpacePoints": ckfSpacePoints,
            "outputSeeds": "gnn-seeds",
            "outputParameters": "gnn-parameters",
            "outputProtoTracks": "gnn-protoTracks-with-parameters",
            "geometry": trackingGeometry,
            "magneticField": field,
            "spacePointSelection": acts.examples.SeedSpacePointSelection.SpreadTriplet,
        }
        parametersAlg = acts.examples.ProtoTracksToParameters(**parametersConfig)
        s.addAlgorithm(parametersAlg)

        ckfConfig = {
            "level": logLevel,
            "inputProtoTracks": parametersAlg.config.outputProtoTracks,
            "inputMeasurements": "measurements",
            "inputInitialTrackParameters": parametersAlg.config.outputParameters,
            "outputTracks": "ckf-tracks",
            # chi2 cut 20, one measurement per surface
            "measurementSelectorCfg": acts.MeasurementSelector.Config(
                [(acts.GeometryIdentifier(), ([], [20.0], [], [1]))]
            ),
            "trackingGeometry": trackingGeometry,
            "magneticField": field,
            "findTracks": acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
                trackingGeometry, field, logLevel
            ),
        }
        ckfAlg = acts.examples.gnn.TrackFindingFromProtoTracksAlgorithm(**ckfConfig)
        s.addAlgorithm(ckfAlg)

        matchAlg = acts.examples.TrackTruthMatcher(
            level=logLevel,
            inputTracks=ckfAlg.config.outputTracks,
            inputParticles="particles",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="ckf_track_particle_matching",
            outputParticleTrackMatching="ckf_particle_track_matching",
            doubleMatching=False,
        )
        s.addAlgorithm(matchAlg)
        s.addWhiteboardAlias(
            "track_particle_matching", matchAlg.config.outputTrackParticleMatching
        )
        s.addWhiteboardAlias(
            "particle_track_matching", matchAlg.config.outputParticleTrackMatching
        )

        addTrackWriters(
            s,
            name="ckf",
            tracks=ckfAlg.config.outputTracks,
            outputDirRoot=outputDir,
            writeSummary=True,
            writeFinderPerformance=True,
            writeFinderNTuple=True,
            logLevel=logLevel,
        )

    s.run()
    return s


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description="Run GNN track finding with module maps on ATLAS data"
    )

    argparser.add_argument(
        "--inputRootDump",
        type=Path,
        required=True,
        help="Path to the input ROOT dump file (ATLAS Athena format)",
    )
    argparser.add_argument(
        "--moduleMapPath",
        type=str,
        required=True,
        help="Path prefix for module map files (without .doublets.root/.triplets.root)",
    )
    argparser.add_argument(
        "--gnnModel",
        type=Path,
        required=True,
        help="Path to the GNN model file (.pt, .onnx, or .engine)",
    )
    argparser.add_argument(
        "--outputDir",
        type=Path,
        default=Path.cwd(),
        help="Output directory for performance files",
    )
    argparser.add_argument(
        "--events",
        type=int,
        default=1,
        help="Number of events to process",
    )
    argparser.add_argument(
        "--bufferEvents",
        type=int,
        default=None,
        help="Number of events to buffer (improves I/O performance)",
    )
    argparser.add_argument(
        "--useEdgeLayerConnector",
        action="store_true",
        default=False,
        help="Use EdgeLayerConnector instead of CudaTrackBuilding for track building",
    )
    argparser.add_argument(
        "--trackBuilder",
        choices=["junctionremoval", "edgelayerconnector", "dwalk"],
        default="junctionremoval",
        help="Track builder to use for graph segmentation",
    )
    argparser.add_argument(
        "--edgeClassifierCut",
        type=float,
        default=None,
        help="Override the edge classifier score cut",
    )
    argparser.add_argument(
        "--pixelOnly",
        action="store_true",
        default=False,
        help="Run on pixel space points only (for pixel-only GNN models)",
    )
    argparser.add_argument(
        "--numFeatures",
        type=int,
        choices=[4, 12],
        default=12,
        help="Number of node features the model expects: 4 (space point R, Phi, Z, "
        "Eta) or 12 (space point + 2 clusters)",
    )
    argparser.add_argument("--debug", action="store_true")

    args = argparser.parse_args()

    runGNN4ITk(
        inputRootDump=args.inputRootDump,
        moduleMapPath=args.moduleMapPath,
        gnnModel=args.gnnModel,
        outputDir=args.outputDir,
        events=args.events,
        bufferEvents=args.bufferEvents,
        useEdgeLayerConnector=args.useEdgeLayerConnector,
        trackBuilder=args.trackBuilder,
        edgeClassifierCut=args.edgeClassifierCut,
        pixelOnly=args.pixelOnly,
        numFeatures=args.numFeatures,
        logLevel=acts.logging.DEBUG if args.debug else acts.logging.INFO,
    )
