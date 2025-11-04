#!/usr/bin/env python3

"""
Example: GNN track finding with metric learning.

Full simulation pipeline using metric learning (embedding + KNN) for graph construction.
Typical workflow for ODD and generic detectors.

All parameters are hardcoded except file paths. Users should copy and modify this script
for their specific needs.
"""

from pathlib import Path
import os
import sys

import acts
import acts.examples
from acts.examples.reconstruction import addGnn
from acts import UnitConstants as u
from digitization import runDigitization


def runGnnMetricLearning(
    trackingGeometry,
    field,
    outputDir,
    digiConfigFile,
    geometrySelection,
    modelDir,
    backend="torch",
    outputRoot=False,
    outputCsv=False,
    s=None,
):
    """
    Run complete simulation + GNN tracking with metric learning.

    This example shows the full chain from particle gun to track finding.
    All GNN parameters are hardcoded for reproducibility.

    Args:
        trackingGeometry: Tracking geometry
        field: Magnetic field
        outputDir: Output directory for ROOT/CSV files
        digiConfigFile: Digitization config file path
        geometrySelection: Geometry selection JSON file path
        modelDir: Directory containing GNN models
        backend: "torch" or "onnx"
        outputRoot: Write ROOT output
        outputCsv: Write CSV output
        s: Optional sequencer (creates new one if None)
    """
    # Run simulation + digitization
    s = runDigitization(
        trackingGeometry,
        field,
        outputDir,
        digiConfigFile=digiConfigFile,
        particlesInput=None,
        outputRoot=outputRoot,
        outputCsv=outputCsv,
        s=s,
    )

    # Configure GNN stages for metric learning workflow
    # All parameters hardcoded based on standard metric learning configuration

    # Stage 1: Graph construction via metric learning
    graphConstructorConfig = {
        "level": acts.logging.INFO,
        "modelPath": str(Path(modelDir) / "embed.pt"),
        "embeddingDim": 8,
        "rVal": 1.6,
        "knnVal": 100,
        "selectedFeatures": [0, 1, 2],  # R, Phi, Z
    }
    graphConstructor = acts.examples.TorchMetricLearning(**graphConstructorConfig)

    # Stage 2: Two-stage edge classification (filter + GNN)
    if backend == "torch":
        filterConfig = {
            "level": acts.logging.INFO,
            "modelPath": str(Path(modelDir) / "filter.pt"),
            "cut": 0.01,
            "nChunks": 5,
            "undirected": False,
            "selectedFeatures": [0, 1, 2],
        }
        gnnConfig = {
            "level": acts.logging.INFO,
            "modelPath": str(Path(modelDir) / "gnn.pt"),
            "cut": 0.5,
            "nChunks": 5,
            "undirected": True,
            "selectedFeatures": [0, 1, 2],
        }
        edgeClassifiers = [
            acts.examples.TorchEdgeClassifier(**filterConfig),
            acts.examples.TorchEdgeClassifier(**gnnConfig),
        ]
    elif backend == "onnx":
        filterConfig = {
            "level": acts.logging.INFO,
            "modelPath": str(Path(modelDir) / "filtering.onnx"),
            "cut": 0.01,
        }
        gnnConfig = {
            "level": acts.logging.INFO,
            "modelPath": str(Path(modelDir) / "gnn.onnx"),
            "cut": 0.5,
        }
        edgeClassifiers = [
            acts.examples.OnnxEdgeClassifier(**filterConfig),
            acts.examples.OnnxEdgeClassifier(**gnnConfig),
        ]
    else:
        raise ValueError(f"Unknown backend: {backend}")

    # Stage 3: CPU track building
    trackBuilderConfig = {
        "level": acts.logging.INFO,
    }
    trackBuilder = acts.examples.BoostTrackBuilding(**trackBuilderConfig)

    # Node features: Standard 3 features (R, Phi, Z)
    nodeFeatures = [
        acts.examples.NodeFeature.R,
        acts.examples.NodeFeature.Phi,
        acts.examples.NodeFeature.Z,
    ]
    featureScales = [1.0, 1.0, 1.0]

    # Add GNN tracking
    addGnn(
        s,
        graphConstructor=graphConstructor,
        edgeClassifiers=edgeClassifiers,
        trackBuilder=trackBuilder,
        nodeFeatures=nodeFeatures,
        featureScales=featureScales,
        trackingGeometry=trackingGeometry,
        geometrySelection=geometrySelection,
        outputDirRoot=outputDir if outputRoot else None,
        logLevel=acts.logging.INFO,
    )

    s.run()
    return s


if __name__ == "__main__":
    # Parse backend from command line
    backend = "torch"
    if "onnx" in sys.argv:
        backend = "onnx"
    elif "torch" in sys.argv:
        backend = "torch"

    # Setup detector (Generic detector for this example)
    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()

    # Magnetic field
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    # Configuration files
    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    geometrySelection = srcdir / "Examples/Configs/generic-seeding-config.json"
    assert geometrySelection.exists(), f"File not found: {geometrySelection}"

    digiConfigFile = srcdir / "Examples/Configs/generic-digi-smearing-config.json"
    assert digiConfigFile.exists(), f"File not found: {digiConfigFile}"

    # Model directory
    if backend == "torch":
        modelDir = Path.cwd() / "torchscript_models"
        assert (modelDir / "embed.pt").exists(), f"Model not found: {modelDir}/embed.pt"
        assert (modelDir / "filter.pt").exists(), f"Model not found: {modelDir}/filter.pt"
        assert (modelDir / "gnn.pt").exists(), f"Model not found: {modelDir}/gnn.pt"
    else:
        modelDir = Path.cwd() / "onnx_models"
        # Note: metric learning still uses embed.pt even with ONNX backend
        assert (Path.cwd() / "torchscript_models" / "embed.pt").exists()
        assert (modelDir / "filtering.onnx").exists()
        assert (modelDir / "gnn.onnx").exists()

    # Sequencer
    s = acts.examples.Sequencer(events=2, numThreads=1)
    s.config.logLevel = acts.logging.INFO

    # Output directory
    outputDir = Path(os.getcwd())

    # Run the workflow
    runGnnMetricLearning(
        trackingGeometry,
        field,
        outputDir,
        digiConfigFile,
        geometrySelection,
        modelDir,
        backend,
        outputRoot=True,
        outputCsv=False,
        s=s,
    )
