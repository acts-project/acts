def test_generic_construction(tmp_path):
    outputDir = str(tmp_path)

    #! [Basic propagation example with GenericDetector]
    import os

    import acts
    from acts import UnitConstants as u
    import acts.examples
    from acts.examples.simulation import (
        addParticleGun,
        EtaConfig,
        ParticleConfig,
        MomentumConfig,
    )
    from acts.examples.root import (
        RootPropagationSummaryWriter,
        RootPropagationStepsWriter,
    )

    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()
    s = acts.examples.Sequencer(events=10)

    rnd = acts.examples.RandomNumbers(seed=42)

    for d in detector.contextDecorators():
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

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))
    stepper = acts.EigenStepper(field)

    propagator = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

    propagationAlgorithm = acts.examples.PropagationAlgorithm(
        propagatorImpl=propagator,
        level=acts.logging.INFO,
        sterileLogger=False,
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

    s.addWriter(
        RootPropagationStepsWriter(
            level=acts.logging.INFO,
            collection="propagation_summary",
            filePath=outputDir + "/propagation_steps.root",
        )
    )

    s.run()
    #! [Basic propagation example with GenericDetector]
