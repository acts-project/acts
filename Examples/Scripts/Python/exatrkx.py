#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union

import acts.examples
import acts
from acts import UnitConstants as u


if "__main__" == __name__:
    import os
    from digitization import configureDigitization
    from acts.examples.reconstruction import addExaTrkx

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

    onnxdir = Path(cwd=os.getcwd()) / "onnx_models"
    assert (
        (onnxdir / "embedding.onnx").exists()
        and (onnxdir / "filtering.onnx").exists()
        and (onnxdir / "gnn.onnx").exists()
    )

    s = acts.examples.Sequencer(events=2, numThreads=-1)
    s.config.logLevel = acts.logging.INFO

    rnd = acts.examples.RandomNumbers()
    outputDir = Path(os.getcwd())

    s = configureDigitization(
        trackingGeometry,
        field,
        outputDir,
        inputParticlePath,
        outputRoot=True,
        outputCsv=True,
        s=s,
    )

    s = addExaTrkx(s, trackingGeometry, geometrySelection, onnxdir, outputDir)

    s.run()
