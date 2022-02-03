#!/usr/bin/env python3
import os

import acts
import acts.examples
import datetime
import argparse

from acts.examples import GenericDetector, AlignedDetector


u = acts.UnitConstants


def runPropagation(
    trackingGeometry,
    field,
    outputDir,
    ntracks,
    outputObj=True,
    outputRoot=True,
    covTransport=False,
    stepper=None,
    s=None,
    decorators=[],
):
    s = s or acts.examples.Sequencer(events=100, numThreads=1)

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    nav = acts.Navigator(trackingGeometry=trackingGeometry)

    stepper = stepper or acts.EigenStepper(field)

    print("We're running with:", type(stepper).__name__)
    prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

    alg = acts.examples.PropagationAlgorithm(
        propagatorImpl=prop,
        level=acts.logging.INFO,
        randomNumberSvc=rnd,
        ntests=ntracks,
        sterileLogger=True,
        covarianceTransport=covTransport,
        recordMaterialInteractions=False,
        propagationStepCollection="propagation-steps",
    )

    s.addAlgorithm(alg)

    # Output
    if outputObj:
        s.addWriter(
            acts.examples.ObjPropagationStepsWriter(
                level=acts.logging.INFO,
                collection="propagation-steps",
                outputDir=outputDir + "/obj",
            )
        )

    if outputRoot:
        s.addWriter(
            acts.examples.RootPropagationStepsWriter(
                level=acts.logging.INFO,
                collection="propagation-steps",
                filePath=outputDir + "/propagation_steps.root",
            )
        )

    return s, alg


if "__main__" == __name__:
    matDeco = None
    # matDeco = acts.IMaterialDecorator.fromFile("material.json")
    # matDeco = acts.IMaterialDecorator.fromFile("material.root")

    ## Generic detector: Default
    (
        detector,
        trackingGeometry,
        contextDecorators,
    ) = GenericDetector.create(mdecorator=matDeco)

    ## Alternative: Aligned detector in a couple of modes
    # detector, trackingGeometry, contextDecorators = AlignedDetector.create(
    #     decoratorLogLevel=acts.logging.INFO,
    #     # These parameters need to be tuned so that GC doesn't break
    #     # with multiple threads
    #     iovSize=10,
    #     flushSize=10,
    #     # External alignment store
    #     mode=AlignedDetector.Config.Mode.External,
    #     # OR: Internal alignment storage
    #     # mode=AlignedDetector.Config.Mode.Internal,
    # )

    ## Alternative: DD4hep detector
    # dd4hepCfg = acts.examples.DD4hepDetector.Config()
    # dd4hepCfg.xmlFileNames = [str(getOpenDataDetectorDirectory()/"xml/OpenDataDetector.xml")]
    # detector = acts.examples.DD4hepDetector()
    # trackingGeometry, contextDecorators = detector.finalize(dd4hepCfg, None)

    ## Magnetic field setup: Default: constant 2T longitudinal field
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    ## Alternative: no B field
    # field = acts.NullBField()

    ## Alternative: Analytical solenoid B field, discretized in an interpolated field map
    # solenoid = acts.SolenoidBField(
    #     radius = 1200*u.mm,
    #     length = 6000*u.mm,
    #     bMagCenter = 2*u.T,
    #     nCoils = 1194
    # )
    # field = acts.solenoidFieldMap(
    #     rlim=(0, 1200*u.mm),
    #     zlim=(-5000*u.mm, 5000*u.mm),
    #     nbins=(50, 50),
    #     field=solenoid
    # )

    p = argparse.ArgumentParser()
    p.add_argument("--stepper", "-s", choices=["atlas", "eigen"])
    p.add_argument("--cov", "-c", choices=["yes", "no"])
    p.add_argument("--events", "-n", type=int, default=1000)

    args = p.parse_args()

    steppers = [acts.AtlasStepper(field), acts.EigenStepper(field)]
    if args.stepper is not None:
        steppers = [
            {"atlas": acts.AtlasStepper(field), "eigen": acts.EigenStepper(field)}[
                args.stepper
            ]
        ]

    covs = [True, False]
    if args.cov is not None:
        covs = [args.cov == "yes"]

    ntracks = 1000

    for stepper in steppers:
        for cov in covs:
            s = acts.examples.Sequencer(
                events=args.events, numThreads=1, logLevel=acts.logging.WARNING
            )

            _, alg = runPropagation(
                trackingGeometry,
                field,
                os.getcwd(),
                ntracks=ntracks,
                outputRoot=False,
                outputObj=False,
                covTransport=cov,
                decorators=contextDecorators,
                stepper=stepper,
                s=s,
            )

            start = datetime.datetime.now()
            s.run()
            end = datetime.datetime.now()

            duration = end - start

            tps = args.events * ntracks / duration.total_seconds()

            print("total steps:", alg.nSteps)
            print(
                f"{type(stepper).__name__} (cov={cov}):",
                f"{duration.total_seconds()}s",
                "/",
                round(tps),
                "tracks per second",
            )
