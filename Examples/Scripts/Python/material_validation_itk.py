#!/usr/bin/env python3
import os
import argparse
from pathlib import Path

from acts.examples import Sequencer, RootMaterialTrackWriter

import acts


def runMaterialValidation(
    trackingGeometry,
    decorators,
    field,
    outputDir,
    outputName="propagation-material",
    dumpPropagationSteps=False,
    s=None,
):
    s = s or Sequencer(events=1000, numThreads=-1)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    nav = acts.Navigator(
        trackingGeometry=trackingGeometry,
        resolveSensitive=True,
        resolveMaterial=True,
        resolvePassive=True,
    )

    # stepper = acts.StraightLineStepper()
    stepper = acts.EigenStepper(field)

    prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

    alg = acts.examples.PropagationAlgorithm(
        propagatorImpl=prop,
        level=acts.logging.INFO,
        randomNumberSvc=acts.examples.RandomNumbers(),
        ntests=1000,
        sterileLogger=False,
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

    if dumpPropagationSteps:
        s.addWriter(
            acts.examples.RootPropagationStepsWriter(
                level=acts.logging.INFO,
                collection=alg.config.propagationStepCollection,
                filePath=outputDir + "/propagation_steps.root",
            )
        )

    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Script to run material validation on ITk geometry"
    )
    p.add_argument(
        "geo_dir",
        help="Input directory containing the ITk standalone geometry. Get in touch if you don't have this.",
    )
    p.add_argument("--material", type=str, default="", help="Material file")

    args = p.parse_args()

    geo_example_dir = Path(args.geo_dir)
    assert geo_example_dir.exists(), "Detector example input directory missing"
    assert os.path.exists(
        args.material
    ), "Invalid file path/name in --material. Please check your input!"

    from acts.examples.itk import buildITkGeometry

    detector, trackingGeometry, decorators = buildITkGeometry(
        geo_example_dir, customMaterialFile=args.material
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    runMaterialValidation(
        trackingGeometry, decorators, field, outputDir=os.getcwd()
    ).run()
