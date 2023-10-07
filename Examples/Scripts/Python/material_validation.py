#!/usr/bin/env python3
import os

from acts.examples import Sequencer, RootMaterialTrackWriter

import acts
from acts import (
    UnitConstants as u,
)
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector


def runMaterialValidation(
    trackingGeometry,
    decorators,
    field,
    outputDir,
    outputName="propagation-material",
    s=None,
):
    s = s or Sequencer(events=1000, numThreads=-1)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    nav = acts.Navigator(trackingGeometry=trackingGeometry)

    stepper = acts.StraightLineStepper()
    # stepper = acts.EigenStepper(field)

    prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

    rnd = acts.examples.RandomNumbers(seed=42)

    alg = acts.examples.PropagationAlgorithm(
        propagatorImpl=prop,
        level=acts.logging.INFO,
        randomNumberSvc=rnd,
        ntests=1000,
        sterileLogger=True,
        propagationStepCollection="propagation-steps",
        recordMaterialInteractions=True,
        d0Sigma=0,
        z0Sigma=0,
    )

    s.addAlgorithm(alg)

    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            collection=alg.config.propagationMaterialCollection,
            filePath=os.path.join(outputDir, (outputName + ".root")),
            storeSurface=True,
            storeVolume=True,
        )
    )

    return s


if "__main__" == __name__:
    matDeco = acts.IMaterialDecorator.fromFile("material-map.json")

    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory(), mdecorator=matDeco
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    runMaterialValidation(
        trackingGeometry, decorators, field, outputDir=os.getcwd()
    ).run()
