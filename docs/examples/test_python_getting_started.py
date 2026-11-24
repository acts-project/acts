def test_first_python_run():
    #! [First Python run]
    import acts
    from acts import UnitConstants as u
    import acts.examples
    from acts.examples.simulation import (
        EtaConfig,
        MomentumConfig,
        ParticleConfig,
        addParticleGun,
    )

    s = acts.examples.Sequencer(events=5, numThreads=1)
    addParticleGun(
        s,
        momentumConfig=MomentumConfig(1 * u.GeV, 10 * u.GeV, transverse=True),
        etaConfig=EtaConfig(-2.0, 2.0),
        particleConfig=ParticleConfig(1, acts.PdgParticle.eMuon),
        rnd=acts.examples.RandomNumbers(seed=42),
        printParticles=True,
    )
    s.run()
    #! [First Python run]
