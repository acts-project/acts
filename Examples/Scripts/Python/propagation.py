#!/usr/bin/env python3

import os
import acts
import acts.examples
from acts.examples import GenericDetector, StructureSelector
from acts.examples.alignment import AlignmentDecorator
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import (
    addParticleGun,
    EtaConfig,
    ParticleConfig,
    MomentumConfig,
)

from acts.examples.root import RootPropagationSummaryWriter, RootPropagationStepsWriter

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
        RootPropagationSummaryWriter(
            level=acts.logging.INFO,
            inputSummaryCollection="propagation_summary",
            filePath=outputDir + "/propagation_summary.root",
        )
    )

    if sterileLogger is False:
        s.addWriter(
            RootPropagationStepsWriter(
                level=acts.logging.INFO,
                collection="propagation_summary",
                filePath=outputDir + "/propagation_steps.root",
            )
        )

    return s


if "__main__" == __name__:
    matDeco = None
    contextDecorators = []
    # matDeco = acts.IMaterialDecorator.fromFile("material.json")
    # matDeco = acts.IMaterialDecorator.fromFile("material.root")

    ## Generic detector: Default
    detector = GenericDetector(materialDecorator=matDeco)
    trackingGeometry = detector.trackingGeometry()

    ## Alternative: Aligned Generic detector
    # detector = AlignedGenericDetector(materialDecorator=matDeco)

    ## Alternative: DD4hep detector
    # detector = getOpenDataDetector()
    # trackingGeometry = detector.trackingGeometry()

    ## Alternative: Misaligned DD4hep detector
    # detector = getOpenDataDetector(misaligned=True)
    # trackingGeometry = detector.trackingGeometry()
    # structureSelector = StructureSelector(trackingGeometry)
    # pixelBarrelID = acts.GeometryIdentifier(volume=17)
    # pixelBarrelTransforms = structureSelector.selectedTransforms(
    #     acts.GeometryContext(), pixelBarrelID
    # )
    # alignDecoConfig = AlignmentDecorator.Config()
    # alignDecoConfig.nominalStore = acts.examples.GeoIdAlignmentStore(
    #     pixelBarrelTransforms
    # )

    # gRot = acts.examples.AlignmentGeneratorGlobalRotation()
    # gRot.axis = acts.Vector3(0.0, 0.0, 1.0)
    # gRot.angle = 0.05

    # lShift = acts.examples.AlignmentGeneratorLocalShift()
    # lShift.axisDirection = acts.AxisDirection.AxisZ
    # lShift.shift = 3.0

    # alignDecoConfig.iovGenerators = [((0, 25), lShift), ((25, 50), gRot)]
    # alignDecoConfig.garbageCollection = True
    # alignDecoConfig.gcInterval = 20

    # alignDeco = AlignmentDecorator(alignDecoConfig, acts.logging.VERBOSE)
    # contextDecorators = [alignDeco]

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
        trackingGeometry=trackingGeometry,
        field=field,
        outputDir=os.getcwd() + "/propagation",
        s=None,
        decorators=contextDecorators,
        sterileLogger=True,
    ).run()
