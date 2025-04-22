#!/usr/bin/env python3

import os
import acts
import acts.examples
from acts.examples import GenericDetector, AlignedDetector
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import (
    addParticleGun,
    EtaConfig,
    ParticleConfig,
    MomentumConfig,
)

u = acts.UnitConstants


def runPropagation(
    trackingGeometry,
    field,
    outputDir,
    tracks: int,
    level: acts.logging.Level,
    s=None,
    decorators=[],
    sterileLogger=True,
):
    s = s or acts.examples.Sequencer(events=1000, numThreads=1)

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        ParticleConfig(num=tracks, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
        EtaConfig(-4.0, 4.0),
        MomentumConfig(1 * u.GeV, 100 * u.GeV, transverse=True),
        rnd=rnd,
    )

    trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
        level=acts.logging.WARNING,
        inputParticles="particles_generated",
        outputTrackParameters="params_particles_generated",
    )
    s.addAlgorithm(trkParamExtractor)

    nav = acts.Navigator(trackingGeometry=trackingGeometry, level=level)

    stepper = acts.EigenStepper(field)
    # stepper = acts.AtlasStepper(field)
    # stepper = acts.StraightLineStepper()

    propagator = acts.examples.ConcretePropagator(
        acts.Propagator(stepper, nav, loglevel=level)
    )

    propagationAlgorithm = acts.examples.PropagationAlgorithm(
        propagatorImpl=propagator,
        level=acts.logging.INFO,
        sterileLogger=sterileLogger,
        inputTrackParameters="params_particles_generated",
        outputSummaryCollection="propagation_summary",
    )
    s.addAlgorithm(propagationAlgorithm)

    s.addWriter(
        acts.examples.RootPropagationSummaryWriter(
            level=acts.logging.INFO,
            inputSummaryCollection="propagation_summary",
            filePath=outputDir + "/propagation_summary.root",
        )
    )

    if sterileLogger:
        s.addWriter(
            acts.examples.RootPropagationStepsWriter(
                level=acts.logging.INFO,
                collection="propagation_summary",
                filePath=outputDir + "/propagation_steps.root",
            )
        )

    return s


if "__main__" == __name__:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--events", "-n", type=int, default=1000)
    p.add_argument("--tracks", "-t", type=int, default=1000)
    p.add_argument("--geo", type=str, default="gen1", choices=["gen1", "gen3"])
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()
    matDeco = None
    # matDeco = acts.IMaterialDecorator.fromFile("material.json")
    # matDeco = acts.IMaterialDecorator.fromFile("material.root")

    if args.verbose:
        level = acts.logging.VERBOSE
    else:
        level = acts.logging.INFO

    ## Generic detector: Default
    detector = GenericDetector(
        gen3=args.geo == "gen3",
        materialDecorator=matDeco,
        logLevel=level,
    )

    ## Alternative: Aligned detector in a couple of modes
    # detector = AlignedDetector(
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
    # detector = getOpenDataDetector()
    # trackingGeometry = detector.trackingGeometry()

    trackingGeometry = detector.trackingGeometry()
    contextDecorators = detector.contextDecorators()

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

    os.makedirs(os.getcwd() + "/propagation", exist_ok=True)
    s = acts.examples.Sequencer(events=args.events, numThreads=1)

    runPropagation(
        trackingGeometry,
        field,
        os.getcwd() + "/propagation",
        decorators=contextDecorators,
        sterileLogger=True,
        tracks=args.tracks,
        level=level,
        s=s,
    ).run()
