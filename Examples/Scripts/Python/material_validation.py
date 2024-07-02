#!/usr/bin/env python3

import os
import argparse

import acts
from acts.examples import Sequencer, RootMaterialTrackWriter
from acts.examples.odd import getOpenDataDetector


def runMaterialValidation(
    nevents,
    ntracks,
    trackingGeometry,
    decorators,
    field,
    outputDir,
    outputName="propagation-material",
    s=None,
):
    # Create a sequencer
    s = s or Sequencer(events=nevents, numThreads=-1)

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
        ntests=ntracks,
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
            inputMaterialTracks=alg.config.propagationMaterialCollection,
            filePath=os.path.join(outputDir, (outputName + ".root")),
            storeSurface=True,
            storeVolume=True,
        )
    )

    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser()

    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to process"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=1000, help="Number of tracks per event"
    )
    p.add_argument(
        "-m", "--map", type=str, default="", help="Input file for the material map"
    )
    p.add_argument("-o", "--output", type=str, default="", help="Output file name")

    args = p.parse_args()

    matDeco = acts.IMaterialDecorator.fromFile(args.map)

    detector, trackingGeometry, decorators = getOpenDataDetector(mdecorator=matDeco)

    field = acts.ConstantBField(acts.Vector3(0, 0, 0 * acts.UnitConstants.T))

    runMaterialValidation(
        args.events,
        args.tracks,
        trackingGeometry,
        decorators,
        field,
        outputDir=os.getcwd(),
        outputName=args.output,
    ).run()
