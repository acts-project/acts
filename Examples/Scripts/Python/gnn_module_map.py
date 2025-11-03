#!/usr/bin/env python3

"""
Example: GNN track finding with module maps.

Demonstrates reading pre-simulated ATLAS data and running GNN with module map
(geometry-based) graph construction. Typical workflow for ATLAS ITk.

All parameters are hardcoded except file paths. Users should copy and modify this script
for their specific needs.
"""

from pathlib import Path
import argparse

import acts
import acts.examples
from acts.examples.reconstruction import addGnn

u = acts.UnitConstants


def runGnnModuleMap(
    inputRootDump: Path,
    moduleMapPath: str,
    gnnModel: Path,
    outputDir: Path,
    events: int = 1,
):
    """
    Run GNN tracking with module maps on ATLAS Athena dumps.

    This example shows reading pre-simulated data and running GNN with
    geometry-based graph construction. All GNN parameters hardcoded.

    Args:
        inputRootDump: Path to input ROOT file (ATLAS Athena dump format)
        moduleMapPath: Path prefix for module map files
                      (will load .doublets.root and .triplets.root)
        gnnModel: Path to trained model (.pt, .onnx, or .engine)
        outputDir: Output directory for performance files
        events: Number of events to process
    """
    # Validate inputs
    assert inputRootDump.exists(), f"Input file not found: {inputRootDump}"
    assert Path(moduleMapPath + ".doublets.root").exists(), \
        f"Module map not found: {moduleMapPath}.doublets.root"
    assert Path(moduleMapPath + ".triplets.root").exists(), \
        f"Module map not found: {moduleMapPath}.triplets.root"
    assert gnnModel.exists(), f"Model file not found: {gnnModel}"

    s = acts.examples.Sequencer(events=events, numThreads=1)

    # Read ATLAS Athena ROOT dump
    s.addReader(
        acts.examples.RootAthenaDumpReader(
            level=acts.logging.INFO,
            treename="GNN4ITk",
            inputfiles=[str(inputRootDump)],
            outputSpacePoints="spacepoints",
            outputClusters="clusters",
            outputMeasurements="measurements",
            outputMeasurementParticlesMap="measurement_particles_map",
            outputParticleMeasurementsMap="particle_measurements_map",
            outputParticles="particles",
            skipOverlapSPsPhi=True,
            skipOverlapSPsEta=False,
            absBoundaryTolerance=0.01 * u.mm,
        )
    )

    # Configure GNN stages for module map workflow
    # All parameters hardcoded based on ITk configuration

    # Stage 1: Graph construction via module map
    graphConstructor = acts.examples.ModuleMapCuda(
        moduleMapPath=moduleMapPath,
        rScale=1000.0,
        phiScale=3.141592654,
        zScale=1000.0,
        etaScale=1.0,
        gpuDevice=0,
        gpuBlocks=512,
        moreParallel=True,
    )

    # Stage 2: Single-stage edge classification (auto-detect backend)
    if gnnModel.suffix == ".pt":
        edgeClassifiers = [
            acts.examples.TorchEdgeClassifier(
                modelPath=str(gnnModel),
                cut=0.5,
                useEdgeFeatures=True,
            )
        ]
    elif gnnModel.suffix == ".onnx":
        edgeClassifiers = [
            acts.examples.OnnxEdgeClassifier(
                modelPath=str(gnnModel),
                cut=0.5,
            )
        ]
    elif gnnModel.suffix == ".engine":
        edgeClassifiers = [
            acts.examples.TensorRTEdgeClassifier(
                modelPath=str(gnnModel),
                cut=0.5,
            )
        ]
    else:
        raise ValueError(f"Unsupported model format: {gnnModel.suffix}")

    # Stage 3: GPU track building
    trackBuilder = acts.examples.CudaTrackBuilding(
        useOneBlockImplementation=False,
        doJunctionRemoval=True,
    )

    # Node features: ITk 12-feature configuration (spacepoint + 2 clusters)
    e = acts.examples.NodeFeature
    nodeFeatures = [
        e.R, e.Phi, e.Z, e.Eta,
        e.Cluster1R, e.Cluster1Phi, e.Cluster1Z, e.Cluster1Eta,
        e.Cluster2R, e.Cluster2Phi, e.Cluster2Z, e.Cluster2Eta,
    ]
    featureScales = [1000.0, 3.141592654, 1000.0, 1.0] * 3

    # Add GNN tracking (spacepoints already created by reader, so no trackingGeometry)
    addGnn(
        s,
        graphConstructor=graphConstructor,
        edgeClassifiers=edgeClassifiers,
        trackBuilder=trackBuilder,
        nodeFeatures=nodeFeatures,
        featureScales=featureScales,
        inputSpacePoints="spacepoints",
        inputClusters="clusters",
        outputDirRoot=str(outputDir),
        logLevel=acts.logging.INFO,
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

    args = argparser.parse_args()

    runGnnModuleMap(
        inputRootDump=args.inputRootDump,
        moduleMapPath=args.moduleMapPath,
        gnnModel=args.gnnModel,
        outputDir=args.outputDir,
        events=args.events,
    )
