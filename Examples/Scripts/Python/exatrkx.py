#!/usr/bin/env python3
from pathlib import Path

import acts.examples
import acts
from acts import UnitConstants as u


if "__main__" == __name__:
    import os
    import sys
    from digitization import runDigitization
    from acts.examples.reconstruction import addExaTrkX, ExaTrkXBackend

    backend = ExaTrkXBackend.Torch

    if "onnx" in sys.argv:
        backend = ExaTrkXBackend.Onnx
    if "torch" in sys.argv:
        backend = ExaTrkXBackend.Torch

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    geometrySelection = (
        srcdir
        / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
    )
    assert geometrySelection.exists()

    digiConfigFile = (
        srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    )
    assert digiConfigFile.exists()

    if backend == ExaTrkXBackend.Torch:
        modelDir = Path.cwd() / "torchscript_models"
        assert (modelDir / "embed.pt").exists()
        assert (modelDir / "filter.pt").exists()
        assert (modelDir / "gnn.pt").exists()
    else:
        modelDir = Path.cwd() / "onnx_models"
        assert (modelDir / "embedding.onnx").exists()
        assert (modelDir / "filtering.onnx").exists()
        assert (modelDir / "gnn.onnx").exists()

    s = acts.examples.Sequencer(events=2, numThreads=1)
    s.config.logLevel = acts.logging.INFO

    rnd = acts.examples.RandomNumbers()
    outputDir = Path(os.getcwd())

    s = runDigitization(
        trackingGeometry,
        field,
        outputDir,
        digiConfigFile=digiConfigFile,
        particlesInput=inputParticlePath,
        outputRoot=True,
        outputCsv=True,
        s=s,
    )

    addExaTrkX(
        s,
        trackingGeometry,
        geometrySelection,
        modelDir,
        outputDir,
        backend=backend,
    )

    s.run()
