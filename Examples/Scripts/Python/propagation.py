#!/usr/bin/env python3

import os
import acts
import acts.examples
from acts.examples import GenericDetector, AlignedGenericDetector
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import (
    addParticleGun,
    EtaConfig,
    ParticleConfig,
    MomentumConfig,
)

u = acts.UnitConstants


def runPropagation(
    trackingGeometry, field, outputDir, s=None, decorators=[], sterileLogger=True
):
    s = s or acts.examples.Sequencer(events=100, numThreads=1)

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        ParticleConfig(num=1000, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
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

    nav = acts.Navigator(trackingGeometry=trackingGeometry)

    stepper = acts.EigenStepper(field)
    # stepper = acts.AtlasStepper(field)
    # stepper = acts.StraightLineStepper()

    propagator = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

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
    matDeco = None
    # matDeco = acts.IMaterialDecorator.fromFile("material.json")
    # matDeco = acts.IMaterialDecorator.fromFile("material.root")

    ## Generic detector: Default
    detector = GenericDetector(materialDecorator=matDeco)

    ## Alternative: Aligned Generic detector
    # detector = AlignedGenericDetector(materialDecorator=matDeco)

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

    runPropagation(
        trackingGeometry,
        field,
        os.getcwd() + "/propagation",
        decorators=contextDecorators,
        sterileLogger=True,
    ).run()
