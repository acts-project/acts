import acts

u = acts.UnitConstants


def test_uniform_vertex_generator():
    vertex_generator = acts.examples.UniformVertexGenerator(
        min=acts.Vector4(
            -1 * u.um,
            -2 * u.um,
            -3 * u.mm,
            -4 * u.ns,
        ),
        max=acts.Vector4(
            1 * u.um,
            2 * u.um,
            3 * u.mm,
            4 * u.ns,
        ),
    )

    assert vertex_generator.min[0] == -1 * u.um
    assert vertex_generator.min[1] == -2 * u.um
    assert vertex_generator.min[2] == -3 * u.mm
    assert vertex_generator.min[3] == -4 * u.ns
    assert vertex_generator.max[0] == 1 * u.um
    assert vertex_generator.max[1] == 2 * u.um
    assert vertex_generator.max[2] == 3 * u.mm
    assert vertex_generator.max[3] == 4 * u.ns


def test_additive_vertex_generator():
    g1 = acts.examples.FixedVertexGenerator(
        fixed=acts.Vector4(1 * u.mm, 2 * u.mm, 3 * u.mm, 0)
    )
    g2 = acts.examples.FixedVertexGenerator(
        fixed=acts.Vector4(10 * u.mm, 20 * u.mm, 30 * u.mm, 0)
    )
    additive = acts.examples.AdditiveVertexGenerator(generators=[g1, g2])

    assert len(additive.generators) == 2


def test_lumi_block_vertex_generator():
    gen = acts.examples.LumiBlockVertexGenerator(
        blockSize=1000,
        stddev=acts.Vector4(0.01 * u.mm, 0.01 * u.mm, 1.0 * u.mm, 0),
    )

    assert gen.blockSize == 1000
    assert gen.stddev[0] == 0.01 * u.mm
    assert gen.stddev[1] == 0.01 * u.mm
    assert gen.stddev[2] == 1.0 * u.mm
    assert gen.stddev[3] == 0


def test_examples_fatras_aliases_present():
    assert hasattr(acts.examples, "SimBarcode")
    assert hasattr(acts.examples, "GenerationProcess")
    assert hasattr(acts.examples, "SimulationOutcome")
    assert hasattr(acts.examples, "SimParticleState")

    assert acts.examples.SimBarcode is acts.fatras.Barcode
    assert acts.examples.GenerationProcess is acts.fatras.GenerationProcess
    assert acts.examples.SimulationOutcome is acts.fatras.SimulationOutcome
    assert acts.examples.SimParticleState is acts.fatras.Particle
