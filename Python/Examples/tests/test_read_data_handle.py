import pytest
import acts
import acts.examples

u = acts.UnitConstants

hepmc3 = pytest.importorskip("acts.examples.hepmc3")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class _Noop(acts.examples.IAlgorithm):
    """Minimal algorithm whose only purpose is to serve as a handle parent."""

    def __init__(self):
        acts.examples.IAlgorithm.__init__(self, "_Noop", acts.logging.WARNING)

    def execute(self, context):
        return acts.examples.ProcessCode.SUCCESS


def _make_sequencer(events=3):
    """Sequencer with a particle-gun → HepMC3-converter chain."""
    rng = acts.examples.RandomNumbers(seed=42)
    s = acts.examples.Sequencer(events=events, numThreads=1, logLevel=acts.logging.INFO)

    evGen = acts.examples.EventGenerator(
        level=acts.logging.WARNING,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=1),
                vertex=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(0, 0, 0, 0),
                    mean=acts.Vector4(0, 0, 0, 0),
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(1 * u.GeV, 10 * u.GeV),
                    eta=(-2, 2),
                    numParticles=4,
                ),
            )
        ],
        outputEvent="particle_gun_event",
        randomNumbers=rng,
    )
    s.addReader(evGen)

    converter = hepmc3.HepMC3InputConverter(
        level=acts.logging.WARNING,
        inputEvent="particle_gun_event",
        outputParticles="particles_generated",
        outputVertices="vertices_generated",
        mergePrimaries=False,
    )
    s.addAlgorithm(converter)

    return s


# ---------------------------------------------------------------------------
# Unit-level tests (no sequencer run needed)
# ---------------------------------------------------------------------------


def test_unregistered_type_raises():
    """Constructing a handle for a type without a whiteboard registration fails."""
    dummy = _Noop()
    with pytest.raises(TypeError, match="not registered"):
        acts.examples.ReadDataHandle(dummy, acts.examples.IAlgorithm, "test")


def test_not_initialized_raises():
    """Calling a handle before initialize() raises RuntimeError."""
    dummy = _Noop()
    handle = acts.examples.ReadDataHandle(
        dummy, acts.examples.SimParticleContainer, "InputParticles"
    )
    wb = acts.examples.WhiteBoard(acts.logging.WARNING)
    with pytest.raises(RuntimeError, match="not initialized"):
        handle(wb)


def test_wrong_key_raises():
    class WrongKeyInspector(acts.examples.IAlgorithm):
        def __init__(self):
            super().__init__(name="WrongKeyInspector", level=acts.logging.INFO)
            self.particles = acts.examples.ReadDataHandle(
                self, acts.examples.SimParticleContainer, "InputParticles"
            )
            self.particles.initialize("wrong_key")

        def execute(self, context):
            return acts.examples.ProcessCode.SUCCESS

    s = _make_sequencer(events=3)
    inspector = WrongKeyInspector()
    with pytest.raises(KeyError, match="does not exist"):
        wb = acts.examples.WhiteBoard(acts.logging.WARNING)
        inspector.particles(wb)

    with pytest.raises(
        RuntimeError,
        match="Sequence configuration error: Missing data handle for key 'wrong_key'",
    ):
        s.addAlgorithm(inspector)


# ---------------------------------------------------------------------------
# Integration test
# ---------------------------------------------------------------------------


def test_read_particles_via_handle():
    """A Python IAlgorithm reads SimParticleContainer from the whiteboard."""

    class ParticleInspector(acts.examples.IAlgorithm):
        def __init__(self):
            super().__init__(name="ParticleInspector", level=acts.logging.INFO)
            self.particles = acts.examples.ReadDataHandle(
                self, acts.examples.SimParticleContainer, "InputParticles"
            )
            self.particles.initialize("particles_generated")
            self.seen_events = 0

        def execute(self, context):
            particles = self.particles(context.eventStore)
            assert isinstance(particles, acts.examples.SimParticleContainer)

            self.logger.info(f"Found {len(particles)} particles")
            self.logger.info("Found {} particles <- this also works", len(particles))

            for particle in particles:
                print(particle)

            self.seen_events += 1
            return acts.examples.ProcessCode.SUCCESS

    s = _make_sequencer(events=3)
    inspector = ParticleInspector()
    s.addAlgorithm(inspector)
    s.run()

    assert inspector.seen_events == 3
