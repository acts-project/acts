#!/usr/bin/env python3

from pathlib import Path
import os
import sys

import acts.examples
import acts
from acts.examples.reconstruction import addGnn, addSpacePointsMaking
from acts import UnitConstants as u

from digitization import runDigitization


def runGnnMetricLearning(
    trackingGeometry,
    field,
    outputDir,
    digiConfigFile,
    geometrySelection,
    embedModelPath,
    filterModelPath,
    gnnModelPath,
    outputRoot=False,
    outputCsv=False,
    s=None,
):
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
    
    addSpacePointsMaking(
        s,
        geoSelectionConfigFile=geometrySelection,
        stripGeoSelectionConfigFile=None,
        trackingGeometry=trackingGeometry,
        logLevel=acts.logging.INFO,
    )

    graphConstructorConfig = {
        "level": acts.logging.INFO,
        "modelPath": str(embedModelPath),
        "embeddingDim": 8,
        "rVal": 1.6,
        "knnVal": 100,
        "selectedFeatures": [0, 1, 2],  # R, Phi, Z
    }
    graphConstructor = acts.examples.TorchMetricLearning(**graphConstructorConfig)

    filterConfig = {
        "level": acts.logging.INFO,
        "modelPath": str(filterModelPath),
        "cut": 0.01,
    }
    gnnConfig = {
        "level": acts.logging.INFO,
        "modelPath": str(gnnModelPath),
        "cut": 0.5,
    }

    edgeClassifiers = []

    if filterModelPath.suffix == ".pt":
        edgeClassifiers.append(
            acts.examples.TorchEdgeClassifier(
                **filterConfig,
                nChunks=5,
                undirected=False,
                selectedFeatures=[0,1,2],
            )
        )
    elif filterModelPath.suffix == ".onnx":
        edgeClassifiers.append(acts.examples.OnnxEdgeClassifier(**filterConfig))
    else:
        raise ValueError(f"Unsupported model format: {filterModelPath.suffix}")

    if gnnModelPath.suffix == ".pt":
        edgeClassifiers.append(
            acts.examples.TorchEdgeClassifier(
                **filterConfig,
                undirected=True,
                selectedFeatures=[0,1,2],
            )
        )
    elif gnnModelPath.suffix == ".onnx":
        edgeClassifiers.append(
            acts.examples.OnnxEdgeClassifier(**gnnConfig),
        )
    else:
        raise ValueError(f"Unsupported model format: {filterModelPath.suffix}")

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
        outputDirRoot=outputDir if outputRoot else None,
        logLevel=acts.logging.INFO,
    )

    s.run()


if "__main__" == __name__:
    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    geometrySelection = srcdir / "Examples/Configs/generic-seeding-config.json"
    assert geometrySelection.exists()

    digiConfigFile = srcdir / "Examples/Configs/generic-digi-smearing-config.json"
    assert digiConfigFile.exists()


    # Model paths from ci_models
    ci_models = srcdir / "ci_models/metric_learning"
    if "onnx" in sys.argv:
        embedModelPath = ci_models / "torchscript_models/embed.pt"
        filterModelPath = ci_models / "torchscript_models/filter.pt"
        gnnModelPath = ci_models / "torchscript_models/gnn.pt"
    elif "torch" in sys.argv:
        embedModelPath = ci_models / "torchscript_models/embed.pt"
        filterModelPath = ci_models / "onnx_models/filtering.onnx"
        gnnModelPath = ci_models / "onnx_models/gnn.onnx"
    else:
        raise ValueError("Please specify backend: 'torch' or 'onnx'")

    s = acts.examples.Sequencer(events=2, numThreads=1)
    s.config.logLevel = acts.logging.INFO

    rnd = acts.examples.RandomNumbers()
    outputDir = Path(os.getcwd())

    runGnnMetricLearning(
        trackingGeometry,
        field,
        outputDir,
        digiConfigFile,
        geometrySelection,
        embedModelPath,
        filterModelPath,
        gnnModelPath,
        outputRoot=True,
        outputCsv=False,
        s=s,
    )
