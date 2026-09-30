def test_custom_algorithm():
    #! [Custom algorithm filtering particles by pT]
    import acts
    from acts import UnitConstants as u
    import acts.examples
    from acts.examples.simulation import (
        addParticleGun,
        MomentumConfig,
        EtaConfig,
        ParticleConfig,
    )

    class HighPtParticleFilter(acts.examples.IAlgorithm):
        """Keeps only particles above a transverse-momentum threshold."""

        def __init__(self, ptMin, level):
            acts.examples.IAlgorithm.__init__(self, "HighPtParticleFilter", level)

            self.ptMin = ptMin

            self.inputParticles = acts.examples.ReadDataHandle(
                self, acts.examples.SimParticleContainer, "InputParticles"
            )
            self.inputParticles.initialize("particles_generated")

            self.outputParticles = acts.examples.WriteDataHandle(
                self, acts.examples.SimParticleContainer, "OutputParticles"
            )
            self.outputParticles.initialize("particles_high_pt")

        def execute(self, context):
            particles = self.inputParticles(context.eventStore)

            kept = acts.examples.SimParticleContainer()
            for particle in particles:
                if particle.transverseMomentum >= self.ptMin:
                    kept.insert(particle)

            self.logger.info(
                "Kept {}/{} particles above {} GeV",
                len(kept),
                len(particles),
                self.ptMin / u.GeV,
            )

            self.outputParticles(context, kept)
            return acts.examples.ProcessCode.SUCCESS

    s = acts.examples.Sequencer(events=5, numThreads=1)
    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        MomentumConfig(0.1 * u.GeV, 10 * u.GeV, transverse=True),
        EtaConfig(-2.0, 2.0),
        ParticleConfig(20, acts.PdgParticle.eMuon, randomizeCharge=True),
        rnd=rnd,
    )

    s.addAlgorithm(HighPtParticleFilter(ptMin=1 * u.GeV, level=acts.logging.INFO))

    s.run()
    #! [Custom algorithm filtering particles by pT]
