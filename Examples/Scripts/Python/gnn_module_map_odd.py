#!/usr/bin/env python3

"""
Example: GNN track finding with module maps on ODD.

Full simulation pipeline using module map (geometry-based) graph construction on the
OpenDataDetector. Demonstrates the complete workflow from event generation to track finding.

All parameters are hardcoded except file paths. Users should copy and modify this script
for their specific needs.
"""

from pathlib import Path

import acts
from acts import UnitConstants as u
from acts.examples import Sequencer
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import (
    addPythia8,
    addFatras,
    addDigitization,
    addGenParticleSelection,
    ParticleSelectorConfig
)
from acts.examples.reconstruction import addGnn


def runGnnModuleMapOdd(
    trackingGeometry,
    field,
    geometrySelection,
    digiConfigFile,
    moduleMapPath,
    gnnModel,
    outputDir,
    events=100,
    s=None,
):
    """
    Run GNN tracking with module maps on ODD.

    This example shows the full simulation chain with module map-based GNN.
    All GNN parameters are hardcoded for reproducibility.

    Args:
        trackingGeometry: Tracking geometry
        field: Magnetic field
        geometrySelection: Geometry selection JSON file path
        digiConfigFile: Digitization config file path
        moduleMapPath: Path prefix for module map files
                      (will load .doublets.root and .triplets.root)
        gnnModel: Path to trained model (.pt, .onnx, or .engine)
        outputDir: Output directory for performance files
        events: Number of events to process
        s: Optional sequencer (creates new one if None)
    """
    # Validate inputs
    assert Path(moduleMapPath + ".doublets.root").exists(), \
        f"Module map not found: {moduleMapPath}.doublets.root"
    assert Path(moduleMapPath + ".triplets.root").exists(), \
        f"Module map not found: {moduleMapPath}.triplets.root"
    assert Path(gnnModel).exists(), f"Model file not found: {gnnModel}"

    s = s or Sequencer(events=events, numThreads=1)

    # Random number generator
    rnd = acts.examples.RandomNumbers(seed=42)

    # Pythia8: ttbar events with pile-up 200 (matching full_chain_odd.py)
    addPythia8(
        s,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=200,
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(
                0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
            ),
        ),
        rnd=rnd,
        logLevel=acts.logging.INFO,
    )

    addGenParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            pt=(150 * u.MeV, None),
        ),
    )


    # FATRAS simulation
    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        enableInteractions=True,
        outputDirRoot=None,
        outputDirCsv=None,
        outputDirObj=None,
        logLevel=acts.logging.INFO,
    )

    # Digitization
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
        logLevel=acts.logging.INFO,
    )

    # Configure GNN stages for module map workflow
    # All parameters hardcoded based on ODD + module map configuration

    # Stage 1: Graph construction via module map
    moduleMapConfig = {
        "level": acts.logging.INFO,
        "moduleMapPath": moduleMapPath,
        "rScale": 1000.0,
        "phiScale": 3.141592654,
        "zScale": 1000.0,
        "etaScale": 1.0,
        "gpuDevice": 0,
        "gpuBlocks": 512,
        "moreParallel": True,
    }
    graphConstructor = acts.examples.ModuleMapCuda(**moduleMapConfig)

    # Stage 2: Single-stage edge classification (auto-detect backend)
    gnnModel = Path(gnnModel)
    edgeClassifierConfig = {
        "level": acts.logging.INFO,
        "modelPath": str(gnnModel),
        "cut": 0.5,
    }

    if gnnModel.suffix == ".pt":
        edgeClassifierConfig["useEdgeFeatures"] = True
        edgeClassifiers = [acts.examples.TorchEdgeClassifier(**edgeClassifierConfig)]
    elif gnnModel.suffix == ".onnx":
        edgeClassifiers = [acts.examples.OnnxEdgeClassifier(**edgeClassifierConfig)]
    elif gnnModel.suffix == ".engine":
        edgeClassifiers = [acts.examples.TensorRTEdgeClassifier(**edgeClassifierConfig)]
    else:
        raise ValueError(f"Unsupported model format: {gnnModel.suffix}")

    # Stage 3: GPU track building
    trackBuilderConfig = {
        "level": acts.logging.INFO,
        "useOneBlockImplementation": False,
        "doJunctionRemoval": True,
    }
    trackBuilder = acts.examples.CudaTrackBuilding(**trackBuilderConfig)

    # Node features: ITk 12-feature configuration
    e = acts.examples.NodeFeature
    nodeFeatures = [
        e.R, e.Phi, e.Z, e.Eta,
        e.Cluster1R, e.Cluster1Phi, e.Cluster1Z, e.Cluster1Eta,
        e.Cluster2R, e.Cluster2Phi, e.Cluster2Z, e.Cluster2Eta,
    ]
    featureScales = [1000.0, 3.141592654, 1000.0, 1.0] * 3

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
        inputClusters="clusters",
        outputDirRoot=str(outputDir),
        logLevel=acts.logging.INFO,
    )

    s.run()
    return s


if __name__ == "__main__":
    # Setup detector (ODD)
    detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    # Magnetic field
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    # Configuration files
    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    geometrySelection = srcdir / "Examples/Configs/odd-seeding-config.json"
    assert geometrySelection.exists(), f"File not found: {geometrySelection}"

    digiConfigFile = srcdir / "Examples/Configs/odd-digi-smearing-config.json"
    assert digiConfigFile.exists(), f"File not found: {digiConfigFile}"

    # Hardcoded model paths
    ci_models_odd = srcdir / "ci_models/odd_module_map"
    moduleMapPath = str(ci_models_odd / "module_map_odd_2k_events.1e-03.float")
    gnnModel = str(ci_models_odd / "gnn_odd_module_map.pt")
    outputDir = Path.cwd()
    events = 100

    # Run the workflow
    runGnnModuleMapOdd(
        trackingGeometry=trackingGeometry,
        field=field,
        geometrySelection=geometrySelection,
        digiConfigFile=digiConfigFile,
        moduleMapPath=moduleMapPath,
        gnnModel=gnnModel,
        outputDir=outputDir,
        events=events,
    )
